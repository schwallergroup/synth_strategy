from typing import Tuple, Dict, List
import copy
from rdkit.Chem import AllChem, rdFMCS
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

from pathlib import Path
root_data = Path(__file__).parent.parent

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


NUCLEOPHILIC_SUBSTITUTION_REACTIONS = [
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Williamson Ether Synthesis",
    "S-alkylation of thiols",
    "Buchwald-Hartwig",
    "N-arylation",
    "Mitsunobu",
    "SN2",
]

HETEROCYCLES_OF_INTEREST = [
    "pyridine", "pyrimidine", "pyrazine", "pyridazine", "imidazole",
    "thiazole", "oxazole", "triazole", "tetrazole", "indole",
    "benzimidazole", "quinoline", "isoquinoline", "furan", "thiophene",
    "pyrrole", "oxadiazole", "thiadiazole", "piperidine", "piperazine",
    "morpholine", "thiomorpholine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage nucleophilic substitution involving a benzylic halide on a heterocyclic core and a nitrogen-containing fragment. The reaction must be one of a specific list of named reactions (e.g., N-alkylation, Buchwald-Hartwig). The heterocyclic core must be one of the rings specified in HETEROCYCLES_OF_INTEREST.
    """
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    benzylic_halide_found = False
    late_stage_coupling = False

    print("Starting analysis for benzylic halide coupling strategy...")

    def dfs_traverse(node, depth=0):
        nonlocal benzylic_halide_found, late_stage_coupling, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"\nAnalyzing reaction at depth {depth}: {rsmi}")

                # Check if this is a nucleophilic substitution reaction
                is_nucleophilic_sub = False

                # Check for specific nucleophilic substitution reactions
                for rxn_type in NUCLEOPHILIC_SUBSTITUTION_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_nucleophilic_sub = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        print(f"Found nucleophilic substitution reaction: {rxn_type}")
                        break

                if is_nucleophilic_sub:
                    print(f"Found nucleophilic substitution reaction at depth {depth}")

                    # Identify benzylic halide reactant and nitrogen nucleophile
                    benzylic_halide_reactant = None
                    nitrogen_nucleophile_reactant = None

                    for reactant in reactants:
                        # Check if it's benzylic (connected to aromatic carbon)
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            # More comprehensive benzylic halide check
                            benzylic_pattern = Chem.MolFromSmarts("c[CD1,CD2,CD3][F,Cl,Br,I]")
                            if mol.HasSubstructMatch(benzylic_pattern):
                                findings_json["atomic_checks"]["functional_groups"].append("benzylic halide")
                                print(f"Found potential benzylic halide: {reactant}")

                                # Check if the core is heterocyclic
                                heterocycle_found = False

                                for ring in HETEROCYCLES_OF_INTEREST:
                                    if checker.check_ring(ring, reactant):
                                        heterocycle_found = True
                                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                                        print(f"Found heterocycle ({ring}) in benzylic halide")
                                        break

                                if heterocycle_found:
                                    benzylic_halide_reactant = reactant

                        # Check for nitrogen nucleophile
                        has_amine = False
                        if checker.check_fg("Primary amine", reactant):
                            has_amine = True
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        if checker.check_fg("Secondary amine", reactant):
                            has_amine = True
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                        if checker.check_fg("Tertiary amine", reactant):
                            has_amine = True
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")
                        if checker.check_fg("Aniline", reactant):
                            has_amine = True
                            findings_json["atomic_checks"]["functional_groups"].append("Aniline")

                        if has_amine:
                            print(f"Found nitrogen nucleophile: {reactant}")
                            nitrogen_nucleophile_reactant = reactant

                    # Verify both components are found and the reaction is at a late stage
                    if benzylic_halide_reactant and nitrogen_nucleophile_reactant:
                        print(
                            f"Found both benzylic halide and nitrogen nucleophile at depth {depth}"
                        )
                        # Add structural constraint for co-occurrence
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "scope": "reaction_step",
                                "description": "A single reaction step must be a named nucleophilic substitution and involve specific reactants: one reactant must be a benzylic halide on a specified heterocycle, and another must be a nitrogen nucleophile.",
                                "targets": [
                                    "NUCLEOPHILIC_SUBSTITUTION_REACTIONS",
                                    "HETEROCYCLES_OF_INTEREST",
                                    "benzylic halide",
                                    "nitrogen_nucleophile"
                                ]
                            }
                        })

                        # Verify the product contains the heterocycle but not the halide
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            # Check if product still has the heterocycle
                            heterocycle_in_product = False
                            for ring in HETEROCYCLES_OF_INTEREST:
                                if checker.check_ring(ring, product):
                                    heterocycle_in_product = True
                                    # No need to add to findings_json here, already added for reactant
                                    print(f"Found heterocycle ({ring}) in product")
                                    break

                            # Check if halide is gone in product
                            halide_in_product = False
                            if checker.check_fg("Primary halide", product):
                                halide_in_product = True
                                findings_json["atomic_checks"]["functional_groups"].append("Primary halide")
                            if checker.check_fg("Secondary halide", product):
                                halide_in_product = True
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary halide")
                            if checker.check_fg("Tertiary halide", product):
                                halide_in_product = True
                                findings_json["atomic_checks"]["functional_groups"].append("Tertiary halide")

                            if halide_in_product:
                                print("Halide still present in product")

                            if heterocycle_in_product and not halide_in_product:
                                benzylic_halide_found = True
                                print("Confirmed benzylic halide coupling transformation")
                                # Add structural constraint for negation
                                findings_json["structural_constraints"].append({
                                    "type": "negation",
                                    "details": {
                                        "scope": "reaction_product",
                                        "description": "The product of the qualifying reaction must not contain a halide functional group.",
                                        "targets": [
                                            "Primary halide",
                                            "Secondary halide",
                                            "Tertiary halide"
                                        ]
                                    }
                                })

                                # Check if this is a late stage reaction (depth <= 3)
                                if depth <= 3:
                                    late_stage_coupling = True
                                    print(
                                        f"SUCCESS: Found late-stage benzylic halide coupling at depth {depth}"
                                    )
                                    # Add structural constraint for positional
                                    findings_json["structural_constraints"].append({
                                        "type": "positional",
                                        "details": {
                                            "target": "benzylic_halide_coupling",
                                            "description": "The qualifying reaction must occur at a late stage, defined as within the first four steps of the synthesis from the final product (depth <= 3).",
                                            "position_type": "depth",
                                            "operator": "<=",
                                            "value": 3
                                        }
                                    })
                                else:
                                    print(
                                        f"Found benzylic halide coupling, but not at late stage (depth {depth} > 3)"
                                    )
                            else:
                                print("Product structure doesn't match expected transformation")
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth
            if node["type"] != "reaction": # This means current node is 'chemical'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    result = benzylic_halide_found and late_stage_coupling
    print(
        f"\nFinal result: benzylic_halide_found={benzylic_halide_found}, late_stage_coupling={late_stage_coupling}"
    )
    return result, findings_json
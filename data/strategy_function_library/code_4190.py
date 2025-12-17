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


HETEROCYCLES_OF_INTEREST = [
    "pyridine", "pyrimidine", "pyrazine", "pyridazine", "triazole",
    "tetrazole", "furan", "thiophene", "pyrrole", "imidazole", "oxazole",
    "thiazole", "isoxazole", "isothiazole", "oxadiazole", "thiadiazole",
    "piperidine", "piperazine", "morpholine", "thiomorpholine", "quinoline",
    "isoquinoline", "indole", "benzimidazole", "benzoxazole",
    "benzothiazole", "purine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage (final reaction step) etherification that couples a primary benzylic alcohol with a heterocycle.
    This function specifically identifies Williamson or Mitsunobu ether synthesis reactions.
    The reactant heterocycle must be one of the specific types defined in the HETEROCYCLES_OF_INTEREST list,
    and the product must contain a benzyl ether moiety and one of the same heterocycles.
    """
    etherification_detected = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Define the structural constraints from the original JSON for easy access
    structural_constraints_definitions = [
        {
            "type": "positional",
            "details": {
                "target": [
                    "Williamson Ether Synthesis",
                    "Mitsunobu aryl ether"
                ],
                "position": "last_stage"
            }
        },
        {
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "benzyl alcohol",
                    "heterocycle"
                ],
                "scope": "reactants",
                "within_reaction": [
                    "Williamson Ether Synthesis",
                    "Mitsunobu aryl ether"
                ]
            }
        },
        {
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "benzyl ether",
                    "heterocycle"
                ],
                "scope": "product",
                "within_reaction": [
                    "Williamson Ether Synthesis",
                    "Mitsunobu aryl ether"
                ]
            }
        }
    ]

    def dfs_traverse(node, depth=0):
        nonlocal etherification_detected, findings_json

        if node["type"] == "reaction" and depth <= 1:  # Late stage (final step)
            # Add positional constraint if this is a late-stage reaction
            if depth == 0 and structural_constraints_definitions[0] not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append(structural_constraints_definitions[0])

            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a specific etherification reaction
                is_williamson = checker.check_reaction("Williamson Ether Synthesis", rsmi)
                is_mitsunobu = checker.check_reaction("Mitsunobu aryl ether", rsmi)

                is_etherification = is_williamson or is_mitsunobu

                if is_williamson and "Williamson Ether Synthesis" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Williamson Ether Synthesis")
                if is_mitsunobu and "Mitsunobu aryl ether" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Mitsunobu aryl ether")

                if is_etherification:
                    # Check for benzyl alcohol in reactants
                    has_benzyl_alcohol = False
                    has_heterocycle_reactant = False
                    benzyl_alcohol_reactant = None

                    for reactant in reactants:
                        # Check for benzyl alcohol (primary alcohol attached to aromatic ring)
                        benzyl_mol = Chem.MolFromSmiles(reactant)
                        if benzyl_mol:
                            # Look for -CH2-OH attached to aromatic carbon
                            benzyl_pattern = Chem.MolFromSmarts("c-[CH2]-[OH]")
                            if benzyl_mol.HasSubstructMatch(benzyl_pattern):
                                has_benzyl_alcohol = True
                                benzyl_alcohol_reactant = reactant
                                if "benzyl alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("benzyl alcohol")
                                # No break here, continue to check for heterocycles in other reactants

                    # Check for heterocycles in other reactants
                    for reactant in reactants:
                        if reactant != benzyl_alcohol_reactant:
                            for ring in HETEROCYCLES_OF_INTEREST:
                                if checker.check_ring(ring, reactant):
                                    has_heterocycle_reactant = True
                                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                                    break
                            if has_heterocycle_reactant:
                                break

                    # Verify the connection in the product
                    if has_benzyl_alcohol and has_heterocycle_reactant:
                        # Add co-occurrence constraint for reactants
                        if structural_constraints_definitions[1] not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(structural_constraints_definitions[1])

                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            # Look for benzyl ether pattern
                            benzyl_ether_pattern = Chem.MolFromSmarts("c-[CH2]-[O]")
                            if product_mol.HasSubstructMatch(benzyl_ether_pattern):
                                if "benzyl ether" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("benzyl ether")

                                # Verify the heterocycle is still present in the product
                                has_heterocycle_in_product = False
                                for ring in HETEROCYCLES_OF_INTEREST:
                                    if checker.check_ring(ring, product):
                                        has_heterocycle_in_product = True
                                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                                        break

                                if has_heterocycle_in_product:
                                    # Add co-occurrence constraint for product
                                    if structural_constraints_definitions[2] not in findings_json["structural_constraints"]:
                                        findings_json["structural_constraints"].append(structural_constraints_definitions[2])
                                    etherification_detected = True

        # Continue traversing
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return etherification_detected, findings_json

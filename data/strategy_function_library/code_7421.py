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
    "purine",
    "pyridine",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "triazole",
    "tetrazole",
    "quinoline",
    "isoquinoline",
    "benzimidazole",
    "benzotriazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving sequential nucleophilic aromatic substitution (SNAr)
    reactions on a halogenated heterocyclic scaffold, with stepwise introduction of amine nucleophiles.
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

    # Track key features
    snar_reactions_count = 0
    has_halogenated_heterocycle = False
    all_nucleophiles_are_amines = True
    heterocycle_maintained = True

    # Track heterocycle transformations
    heterocycle_transformations = []

    def dfs_traverse(node, depth=0):
        nonlocal snar_reactions_count, has_halogenated_heterocycle, all_nucleophiles_are_amines, heterocycle_maintained, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            has_heterocycle = False
            heterocycle_type = None
            for ring in HETEROCYCLES_OF_INTEREST:
                if checker.check_ring(ring, mol_smiles):
                    has_heterocycle = True
                    heterocycle_type = ring
                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                    break

            if has_heterocycle:
                has_halide = False
                halide_types = []
                if checker.check_fg("Aromatic halide", mol_smiles):
                    has_halide = True
                    halide_types.append("Aromatic halide")
                if checker.check_fg("Primary halide", mol_smiles):
                    has_halide = True
                    halide_types.append("Primary halide")
                if checker.check_fg("Secondary halide", mol_smiles):
                    has_halide = True
                    halide_types.append("Secondary halide")
                if checker.check_fg("Tertiary halide", mol_smiles):
                    has_halide = True
                    halide_types.append("Tertiary halide")

                for fg in halide_types:
                    if fg not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(fg)

                if has_halide:
                    has_halogenated_heterocycle = True
                    halogen_count = 0
                    mol = Chem.MolFromSmiles(mol_smiles)
                    if mol:
                        for atom in mol.GetAtoms():
                            if atom.GetSymbol() in ["F", "Cl", "Br", "I"]:
                                halogen_count += 1

                    heterocycle_transformations.append(
                        {
                            "depth": depth,
                            "smiles": mol_smiles,
                            "heterocycle": heterocycle_type,
                            "halogen_count": halogen_count,
                        }
                    )

        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                snar_reaction_names = [
                    "heteroaromatic_nuc_sub",
                    "nucl_sub_aromatic_ortho_nitro",
                    "nucl_sub_aromatic_para_nitro",
                    "N-arylation",
                    "Buchwald-Hartwig",
                    "N-arylation_heterocycles",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine"
                ]

                is_snar = False
                detected_snar_name = None
                for r_name in snar_reaction_names:
                    if checker.check_reaction(r_name, rsmi):
                        is_snar = True
                        detected_snar_name = r_name
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        break

                if is_snar:
                    product_has_heterocycle = any(
                        checker.check_ring(ring, product_smiles) for ring in HETEROCYCLES_OF_INTEREST
                    )

                    has_amine = False
                    amine_types = []
                    for reactant in reactants_smiles:
                        if checker.check_fg("Primary amine", reactant):
                            has_amine = True
                            amine_types.append("Primary amine")
                        if checker.check_fg("Secondary amine", reactant):
                            has_amine = True
                            amine_types.append("Secondary amine")
                        if checker.check_fg("Aniline", reactant):
                            has_amine = True
                            amine_types.append("Aniline")
                    
                    for fg in amine_types:
                        if fg not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(fg)

                    if product_has_heterocycle and has_amine:
                        snar_reactions_count += 1
                    elif not has_amine:
                        all_nucleophiles_are_amines = False
                        if {"type": "negation", "details": {"target": "SNAr reaction without an amine nucleophile", "description": "The route is invalid if any SNAr-type reaction occurs where the nucleophile is not an amine."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "SNAr reaction without an amine nucleophile", "description": "The route is invalid if any SNAr-type reaction occurs where the nucleophile is not an amine."}})
                    elif not product_has_heterocycle:
                        heterocycle_maintained = False
                        if {"type": "negation", "details": {"target": "heterocycle_destruction", "description": "The route is invalid if a reaction consumes a heterocycle of interest without it being present in the product."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "heterocycle_destruction", "description": "The route is invalid if a reaction consumes a heterocycle of interest without it being present in the product."}})

                reactant_has_heterocycle = False
                for reactant in reactants_smiles:
                    if any(checker.check_ring(ring, reactant) for ring in HETEROCYCLES_OF_INTEREST):
                        reactant_has_heterocycle = True
                        break

                product_has_heterocycle_in_step = any(
                    checker.check_ring(ring, product_smiles) for ring in HETEROCYCLES_OF_INTEREST
                )

                if reactant_has_heterocycle and not product_has_heterocycle_in_step:
                    heterocycle_maintained = False
                    if {"type": "negation", "details": {"target": "heterocycle_destruction", "description": "The route is invalid if a reaction consumes a heterocycle of interest without it being present in the product."}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "heterocycle_destruction", "description": "The route is invalid if a reaction consumes a heterocycle of interest without it being present in the product."}})

            except Exception as e:
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    sequential_snar = False
    if len(heterocycle_transformations) >= 2:
        heterocycle_transformations.sort(key=lambda x: x["depth"])

        for i in range(len(heterocycle_transformations) - 1):
            current = heterocycle_transformations[i]
            next_mol = heterocycle_transformations[i + 1]

            if next_mol["halogen_count"] > current["halogen_count"]:
                sequential_snar = True
                if {"type": "sequence", "details": {"description": "At least one reaction step involving a heterocycle of interest must show a reduction in halogen count from its precursor to its product (i.e., precursor_halogens > product_halogens)."}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "sequence", "details": {"description": "At least one reaction step involving a heterocycle of interest must show a reduction in halogen count from its precursor to its product (i.e., precursor_halogens > product_halogens)."}})
                break

    strategy_present = (
        snar_reactions_count >= 2
        and has_halogenated_heterocycle
        and all_nucleophiles_are_amines
        and heterocycle_maintained
        and sequential_snar
    )

    if snar_reactions_count >= 2:
        if {"type": "count", "details": {"target": "SNAr reaction with an amine nucleophile on a heterocycle of interest", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "SNAr reaction with an amine nucleophile on a heterocycle of interest", "operator": ">=", "value": 2}})
    
    if has_halogenated_heterocycle:
        if {"type": "co-occurrence", "details": {"targets": ["any_heterocycle_of_interest", "any_halide_functional_group"], "description": "At least one molecule in the route must be a halogenated heterocycle from the lists of interest."}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["any_heterocycle_of_interest", "any_halide_functional_group"], "description": "At least one molecule in the route must be a halogenated heterocycle from the lists of interest."}})

    return strategy_present, findings_json

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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving sequential formation and transformation
    of small (3-membered) heterocyclic rings like aziridines, epoxides, and thiiranes.
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

    # Track the presence of small rings and their transformations
    aziridine_formed = False
    epoxide_formed = False
    thiirane_formed = False
    ring_transformations = 0
    stereocenter_inversions = 0

    # Track the sequence of transformations
    transformation_sequence = []

    def dfs_traverse(node, depth=0):
        nonlocal aziridine_formed, epoxide_formed, thiirane_formed, ring_transformations, stereocenter_inversions, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0]
                products_smiles = rsmi.split(">")[-1]

                try:
                    # Check for aziridine formation/transformation
                    reactant_aziridines = checker.check_ring("aziridine", reactants_smiles)
                    product_aziridines = checker.check_ring("aziridine", products_smiles)

                    if product_aziridines and not reactant_aziridines:
                        print(f"Aziridine formation detected at depth {depth}")
                        aziridine_formed = True
                        transformation_sequence.append(("aziridine_formation", depth))
                        if "aziridine" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("aziridine")
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                    elif reactant_aziridines and not product_aziridines:
                        print(f"Aziridine opening detected at depth {depth}")
                        ring_transformations += 1
                        transformation_sequence.append(("aziridine_opening", depth))
                        if "aziridine" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("aziridine")
                        if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")

                    # Check for epoxide formation/transformation
                    reactant_epoxides = checker.check_ring("oxirane", reactants_smiles)
                    product_epoxides = checker.check_ring("oxirane", products_smiles)

                    if product_epoxides and not reactant_epoxides:
                        print(f"Epoxide formation detected at depth {depth}")
                        epoxide_formed = True
                        transformation_sequence.append(("epoxide_formation", depth))
                        if "oxirane" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("oxirane")
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                    elif reactant_epoxides and not product_epoxides:
                        print(f"Epoxide opening detected at depth {depth}")
                        ring_transformations += 1
                        transformation_sequence.append(("epoxide_opening", depth))
                        if "oxirane" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("oxirane")
                        if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")

                    # Check for thiirane formation/transformation
                    reactant_thiiranes = checker.check_ring("thiirane", reactants_smiles)
                    product_thiiranes = checker.check_ring("thiirane", products_smiles)

                    if product_thiiranes and not reactant_thiiranes:
                        print(f"Thiirane formation detected at depth {depth}")
                        thiirane_formed = True
                        transformation_sequence.append(("thiirane_formation", depth))
                        if "thiirane" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("thiirane")
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                    elif reactant_thiiranes and not product_thiiranes:
                        print(f"Thiirane opening detected at depth {depth}")
                        ring_transformations += 1
                        transformation_sequence.append(("thiirane_opening", depth))
                        if "thiirane" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("thiirane")
                        if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")

                    # Check for reactions that typically involve stereochemical inversion
                    if checker.check_reaction("Mitsunobu esterification", rsmi):
                        print(f"Mitsunobu esterification detected at depth {depth}")
                        stereocenter_inversions += 1
                        transformation_sequence.append(("stereocenter_inversion", depth))
                        if "Mitsunobu esterification" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Mitsunobu esterification")
                    
                    if checker.check_reaction("Mitsunobu aryl ether", rsmi):
                        print(f"Mitsunobu aryl ether detected at depth {depth}")
                        stereocenter_inversions += 1
                        transformation_sequence.append(("stereocenter_inversion", depth))
                        if "Mitsunobu aryl ether" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Mitsunobu aryl ether")

                    if checker.check_reaction("Ring opening of epoxide with amine", rsmi):
                        print(f"Ring opening of epoxide with amine detected at depth {depth}")
                        stereocenter_inversions += 1
                        transformation_sequence.append(("stereocenter_inversion", depth))
                        if "Ring opening of epoxide with amine" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("Ring opening of epoxide with amine")

                except Exception as e:
                    print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Analyze the transformation sequence
    has_formation_opening_pair = False

    # Check for any formation followed by opening (in any order)
    formations = [t for t in transformation_sequence if "formation" in t[0]]
    openings = [t for t in transformation_sequence if "opening" in t[0]]

    if formations and openings:
        has_formation_opening_pair = True
        # Add structural constraint for co-occurrence if detected
        if {"type": "co-occurrence", "details": {"targets": ["ring_formation", "ring_destruction"]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["ring_formation", "ring_destruction"]}})

    # Determine if the strategy is present based on observed patterns
    strategy_present = (
        (aziridine_formed or epoxide_formed or thiirane_formed)
        and ring_transformations >= 1
        and stereocenter_inversions >= 1
        and has_formation_opening_pair
    )

    # Add structural constraint for stereoinvertive_reaction count if detected
    if stereocenter_inversions >= 1:
        if {"type": "count", "details": {"target": "stereoinvertive_reaction", "operator": ">=", "value": 1}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "stereoinvertive_reaction", "operator": ">=", "value": 1}})

    print(f"Small ring interconversion strategy detected: {strategy_present}")
    print(f"Aziridine formed: {aziridine_formed}")
    print(f"Epoxide formed: {epoxide_formed}")
    print(f"Thiirane formed: {thiirane_formed}")
    print(f"Ring transformations: {ring_transformations}")
    print(f"Stereocenter inversions: {stereocenter_inversions}")
    print(f"Transformation sequence: {transformation_sequence}")
    print(f"Has formation-opening pair: {has_formation_opening_pair}")

    return strategy_present, findings_json
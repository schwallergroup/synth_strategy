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
    This function detects if the synthetic route involves late-stage functionalization
    with an acetylaminoethoxy group.
    """
    print("Starting late_stage_functionalization_strategy analysis")
    late_stage_functionalization_detected = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_functionalization_detected, findings_json

        # Print current node type and depth for debugging
        print(f"Traversing node type: {node['type']}, depth: {depth}")

        if node["type"] == "reaction":
            # Use provided depth or calculate from traversal
            node_depth = depth
            if "metadata" in node and "depth" in node["metadata"]:
                node_depth = node["metadata"]["depth"]

            print(f"Reaction node at depth {node_depth}")

            # Check reactions in the late stage (depth 0, 1, 2, or 3)
            if node_depth <= 3:
                # Add positional constraint if met
                if {"type": "positional", "details": {"target": "acetylaminoethoxy_formation", "position": "depth <= 3"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "acetylaminoethoxy_formation", "position": "depth <= 3"}})

                try:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    print(f"Analyzing reaction at depth {node_depth}")
                    print(f"Product: {product_smiles}")
                    print(f"Reactants: {', '.join(reactants_smiles)}")

                    # Check if product has both ether and secondary amide groups
                    has_ether_product = checker.check_fg("Ether", product_smiles)
                    if has_ether_product:
                        if "Ether" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ether")

                    has_secondary_amide_product = checker.check_fg(
                        "Secondary amide", product_smiles
                    )
                    if has_secondary_amide_product:
                        if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")

                    # Check if any reactant has both ether and secondary amide groups
                    has_both_in_reactants = False
                    has_acetylaminoethanol_reactant = False

                    for reactant in reactants_smiles:
                        # Check if both functional groups exist in any reactant
                        if checker.check_fg("Ether", reactant) and checker.check_fg(
                            "Secondary amide", reactant
                        ):
                            has_both_in_reactants = True
                            # This is a negative condition, so we don't add it to findings_json directly here

                        # Check for 2-acetylaminoethanol pattern in reactants
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            acetylaminoethanol_pattern = Chem.MolFromSmarts("CC(=O)NCCO")
                            if reactant_mol.HasSubstructMatch(acetylaminoethanol_pattern):
                                print("Found 2-acetylaminoethanol reactant")
                                has_acetylaminoethanol_reactant = True
                                if "2-acetylaminoethanol" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("2-acetylaminoethanol")

                    print(f"Product has ether: {has_ether_product}")
                    print(f"Product has secondary amide: {has_secondary_amide_product}")
                    print(f"Any reactant has both groups: {has_both_in_reactants}")

                    # Check for acetylaminoethoxy-like structure in product
                    if has_ether_product and has_secondary_amide_product:
                        # Convert to RDKit molecules for more detailed analysis
                        product_mol = Chem.MolFromSmiles(product_smiles)

                        if product_mol:
                            # Check for specific acetylaminoethoxy pattern
                            acetylaminoethoxy_pattern = Chem.MolFromSmarts("CC(=O)NCC[O]")

                            if product_mol.HasSubstructMatch(acetylaminoethoxy_pattern):
                                print("Product has acetylaminoethoxy group")
                                if "acetylaminoethoxy" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("acetylaminoethoxy")

                                # If we found the pattern in product and it was not present in reactants
                                # (approximated by checking for both FGs) or if the specific precursor
                                # 2-acetylaminoethanol was used, it's a valid formation event.
                                if (not has_both_in_reactants) or has_acetylaminoethanol_reactant:
                                    print(
                                        "Confirmed late-stage functionalization with acetylaminoethoxy group"
                                    )
                                    late_stage_functionalization_detected = True
                                    if "acetylaminoethoxy_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                        findings_json["atomic_checks"]["named_reactions"].append("acetylaminoethoxy_formation")

                                    # Add the conditional structural constraint if met
                                    conditional_constraint = {
                                        "type": "conditional",
                                        "details": {
                                            "description": "The formation event is valid if one of the two following precursor conditions is met.",
                                            "conditions": [
                                                {
                                                    "type": "negation",
                                                    "details": {
                                                        "target": "co-occurrence of Ether and Secondary amide",
                                                        "scope": "single_reactant"
                                                    }
                                                },
                                                {
                                                    "type": "co-occurrence",
                                                    "details": {
                                                        "targets": [
                                                            "2-acetylaminoethanol"
                                                        ],
                                                        "scope": "reactants"
                                                    }
                                                }
                                            ],
                                            "logic": "OR"
                                        }
                                    }
                                    if conditional_constraint not in findings_json["structural_constraints"]:
                                        findings_json["structural_constraints"].append(conditional_constraint)

                except Exception as e:
                    print(f"Error processing reaction: {e}")
                    traceback.print_exc() # Print full traceback for debugging

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only from chemical to reaction node
            # Depth remains the same from reaction to chemical node
            new_depth = depth
            if node['type'] != 'reaction': # If current node is chemical, depth increases for reaction child
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Final result: {late_stage_functionalization_detected}")
    return late_stage_functionalization_detected, findings_json

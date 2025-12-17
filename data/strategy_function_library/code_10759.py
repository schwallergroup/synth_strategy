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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the isoxazole structural motif is maintained
    throughout the synthesis, focusing on the target molecule and late-stage
    intermediates.
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

    # Track if we found the motif in the target and late-stage intermediates
    target_and_late_stage_have_motifs = True

    # Track all nodes by their SMILES for depth calculation
    nodes_by_smiles = {}

    # Track products of reactions to identify main synthetic pathway
    reaction_products = set()

    # First pass: collect all nodes and identify reaction products
    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            nodes_by_smiles[node["smiles"]] = depth
        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # Extract product from reaction SMILES
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                product = rsmi.split(">")[-1]
                reaction_products.add(product)
                print(f"Identified reaction product: {product}")
            except Exception as e:
                print(f"Error extracting product from reaction: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Collect all nodes with their depths
    dfs_traverse(route)

    # Second pass: check motifs in target and late-stage intermediates
    def check_motifs(node):
        nonlocal target_and_late_stage_have_motifs, findings_json

        if node["type"] == "mol" and not node.get("in_stock", False):
            depth = nodes_by_smiles.get(node["smiles"], 999)

            # Only check target (depth=0) and late-stage intermediates (depth<=2)
            # that are products of reactions (part of main synthetic pathway)
            if depth <= 2 and (depth == 0 or node["smiles"] in reaction_products):
                print(f"Checking molecule at depth {depth}: {node['smiles']}")

                # Check for isoxazole ring
                has_isoxazole = checker.check_ring("isoxazole", node["smiles"])

                if has_isoxazole:
                    if "isoxazole" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("isoxazole")
                else:
                    target_and_late_stage_have_motifs = False
                    print(
                        f"Molecule at depth {depth} missing the isoxazole motif: {node['smiles']}"
                    )
                    print(f"  Has isoxazole: {has_isoxazole}")

        # Process children
        for child in node.get("children", []):
            check_motifs(child)

    # Start traversal to check motifs
    check_motifs(route)

    if target_and_late_stage_have_motifs:
        print(
            "Target molecule and late-stage intermediates contain the required motif(s)"
        )
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "isoxazole",
                "position": "late_stage",
                "condition": "must_be_present_in_all_molecules_at_depth_<=_2"
            }
        })
    else:
        print("One or more key molecules are missing the required motif(s)")

    return target_and_late_stage_have_motifs, findings_json
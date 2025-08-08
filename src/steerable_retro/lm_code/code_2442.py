#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
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

root_data = "/home/dparm/steerable_retro/data"

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


def main(route):
    """
    Detects if the route contains a late-stage amide coupling (at depth 0 or 1).
    """
    found_late_amide = False

    def calculate_depth(node, root):
        """Calculate depth of a node if not already provided"""
        if "depth" in node:
            return node["depth"]

        # Simple depth calculation based on distance from root
        depth = 0
        current = node
        queue = [(root, 0)]

        while queue:
            current_node, current_depth = queue.pop(0)

            if current_node == current:
                depth = current_depth
                break

            for child in current_node.get("children", []):
                queue.append((child, current_depth + 1))

        return depth

    def dfs_traverse(node, current_depth=0):
        nonlocal found_late_amide

        # Store depth for future reference
        if "depth" not in node:
            node["depth"] = current_depth

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            depth = node.get("depth", current_depth)

            if depth <= 1:  # Late stage (depth 0 or 1)
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for amide coupling reactions using the checker
                amide_coupling_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                ]

                for reaction_type in amide_coupling_reactions:
                    try:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(
                                f"Found late-stage amide coupling: {reaction_type} at depth {depth}"
                            )
                            found_late_amide = True
                            break
                    except Exception as e:
                        print(f"Error checking reaction {reaction_type}: {e}")

                # If no specific reaction type matched, check for amide formation more generally
                if not found_late_amide:
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        has_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants)
                        has_amine = any(
                            checker.check_fg("Primary amine", r) for r in reactants
                        ) or any(checker.check_fg("Secondary amine", r) for r in reactants)
                        forms_amide = (
                            checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                        )

                        if has_acid and has_amine and forms_amide:
                            print(f"Found late-stage amide coupling (generic) at depth {depth}")
                            found_late_amide = True
                    except Exception as e:
                        print(f"Error in generic amide coupling check: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)
    return found_late_amide

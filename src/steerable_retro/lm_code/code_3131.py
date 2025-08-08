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
    This function detects if chlorophenyl groups are maintained throughout the synthesis.
    It checks if the chlorophenyl group, once introduced, is preserved in all subsequent steps.

    Returns:
        bool: True if chlorophenyl groups are maintained throughout synthesis,
              False if they are lost at any point.
    """

    def dfs_traverse(node, depth=0):
        # For molecule nodes, check if it contains a chlorophenyl group
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_chlorophenyl = checker.check_fg("Aromatic halide", mol_smiles)

            # If this is a final product (no children), it must have chlorophenyl
            if not node.get("children") and not has_chlorophenyl:
                print(f"Final product without chlorophenyl found: {mol_smiles}")
                return False

            # For intermediate molecules with children, check if chlorophenyl is maintained
            if has_chlorophenyl:
                # If molecule has chlorophenyl, all its children reactions must preserve it
                for child in node.get("children", []):
                    if not dfs_traverse(child, depth + 1):
                        return False
                return True
            else:
                # If molecule doesn't have chlorophenyl, continue checking children
                for child in node.get("children", []):
                    if dfs_traverse(child, depth + 1):
                        return True
                return False

        # For reaction nodes, check if chlorophenyl is maintained from product to reactants
        elif node["type"] == "reaction":
            try:
                if "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    product = rsmi.split(">")[-1]

                    # If product has chlorophenyl, at least one reactant must have it too
                    if checker.check_fg("Aromatic halide", product):
                        # Check if any child (reactant) has chlorophenyl
                        chlorophenyl_maintained = False
                        for child in node.get("children", []):
                            if dfs_traverse(child, depth + 1):
                                chlorophenyl_maintained = True
                                break

                        if not chlorophenyl_maintained:
                            print(f"Chlorophenyl appeared without source in reaction: {rsmi}")
                            return False
                        return True
                    else:
                        # If product doesn't have chlorophenyl, no need to check reactants
                        return False
            except Exception as e:
                print(f"Error processing reaction: {e}")
                return False

            # Process children for reactions without rsmi
            for child in node.get("children", []):
                if dfs_traverse(child, depth + 1):
                    return True
            return False

        # Default case (shouldn't reach here)
        return False

    # Start traversal from the root (target molecule)
    return dfs_traverse(route)

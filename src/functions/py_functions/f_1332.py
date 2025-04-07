#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects a synthesis involving a fluoroalkoxy substituent
    (O-CH2-CH2-CH2-F) that remains unchanged throughout the synthesis
    """
    # Track if we've found a fluoroalkoxy group that persists through the synthesis
    fluoroalkoxy_preserved = False

    # Define a custom function to check for the fluoroalkoxy group
    # since it's not in the standard functional groups list
    def has_fluoroalkoxy(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False
        pattern = Chem.MolFromSmarts("[O][CH2][CH2][CH2][F]")
        return mol.HasSubstructMatch(pattern)

    # Find the final product (root node)
    if route["type"] == "mol" and has_fluoroalkoxy(route["smiles"]):
        print(f"Final product contains fluoroalkoxy group: {route['smiles']}")

        # Now check if this group is preserved throughout the synthesis
        def trace_fluoroalkoxy(node, depth=0):
            nonlocal fluoroalkoxy_preserved

            if node["type"] == "mol":
                mol_has_group = has_fluoroalkoxy(node["smiles"])
                indent = "  " * depth
                print(
                    f"{indent}Checking molecule: {node['smiles']} - Has fluoroalkoxy: {mol_has_group}"
                )

                # If we reach a starting material with the group, it's preserved
                if mol_has_group and node.get("in_stock", False):
                    print(f"{indent}Found starting material with fluoroalkoxy group")
                    fluoroalkoxy_preserved = True
                    return True

                # If molecule doesn't have the group but isn't a starting material,
                # check its precursors
                if not mol_has_group and not node.get("in_stock", False):
                    return False

            # For reaction nodes or molecules with the group that aren't starting materials
            if "children" in node:
                all_children_have_group = True
                for child in node["children"]:
                    # For reaction nodes, we need at least one child to have the group
                    # For molecule nodes, all children must preserve the group
                    child_result = trace_fluoroalkoxy(child, depth + 1)
                    if node["type"] == "reaction":
                        if child_result:
                            return True
                    else:  # molecule node
                        all_children_have_group = (
                            all_children_have_group and child_result
                        )

                if node["type"] == "mol":
                    return all_children_have_group

            return False

        # Start tracing from the root
        trace_fluoroalkoxy(route)
    else:
        print(f"Final product does not contain fluoroalkoxy group")

    return fluoroalkoxy_preserved

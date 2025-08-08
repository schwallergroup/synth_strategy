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


def main(route):
    """
    This function detects late-stage introduction of trifluoromethyl group.
    """
    cf3_introduced = False
    introduction_depth = -1

    def dfs_traverse(node):
        nonlocal cf3_introduced, introduction_depth

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Extract depth from ID
            depth_match = re.search(r"Depth: (\d+)", node.get("metadata", {}).get("ID", ""))
            current_depth = int(depth_match.group(1)) if depth_match else -1

            # Check for CF3 group
            cf3_pattern = Chem.MolFromSmarts("[C]([F])([F])[F]")
            product_mol = Chem.MolFromSmiles(product)

            if product_mol and product_mol.HasSubstructMatch(cf3_pattern):
                # Check if pattern is not in all reactants (meaning it was introduced)
                cf3_in_all_reactants = True
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and not reactant_mol.HasSubstructMatch(cf3_pattern):
                        cf3_in_all_reactants = False
                        break

                if not cf3_in_all_reactants:
                    cf3_introduced = True
                    introduction_depth = current_depth
                    print(f"Detected CF3 introduction at depth {current_depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Consider it late introduction if it happens in the first half of the synthesis
    # (remember depth 0 is the final step)
    return cf3_introduced and introduction_depth <= 2

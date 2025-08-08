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
    This function detects if a trifluoromethyl group is introduced in the late stage of synthesis.
    """
    trifluoromethyl_introduced = False

    def dfs_traverse(node, depth=0):
        nonlocal trifluoromethyl_introduced

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product has trifluoromethyl
            product_mol = Chem.MolFromSmiles(product_smiles)
            trifluoromethyl_pattern = Chem.MolFromSmarts("[C]([F])([F])[F]")

            if (
                product_mol
                and trifluoromethyl_pattern
                and product_mol.HasSubstructMatch(trifluoromethyl_pattern)
            ):
                # Check if any reactant doesn't have trifluoromethyl
                reactant_without_cf3 = False
                for r_smiles in reactants_smiles:
                    r_mol = Chem.MolFromSmiles(r_smiles)
                    if r_mol and not r_mol.HasSubstructMatch(trifluoromethyl_pattern):
                        reactant_without_cf3 = True
                        break

                if reactant_without_cf3:
                    print(f"Late-stage trifluoromethyl introduction detected in reaction: {rsmi}")
                    trifluoromethyl_introduced = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return trifluoromethyl_introduced

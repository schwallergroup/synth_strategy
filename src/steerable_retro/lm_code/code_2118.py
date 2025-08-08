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
    This function detects if the synthetic route involves manipulation of halogen atoms,
    specifically looking for introduction or removal of iodine while maintaining fluorine atoms.
    """
    aryl_iodide_pattern = Chem.MolFromSmarts("c[I]")
    has_halogen_manipulation = False

    def dfs_traverse(node):
        nonlocal has_halogen_manipulation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Check for iodine addition or removal
                    if reactants_mol and product_mol:
                        reactant_has_iodine = reactants_mol.HasSubstructMatch(aryl_iodide_pattern)
                        product_has_iodine = product_mol.HasSubstructMatch(aryl_iodide_pattern)

                        if reactant_has_iodine != product_has_iodine:
                            print(f"Detected halogen manipulation in reaction: {rsmi}")
                            has_halogen_manipulation = True
                except:
                    print(f"Error processing reaction SMILES: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_halogen_manipulation

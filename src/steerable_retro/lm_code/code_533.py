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
    This function detects if the synthesis includes a nitro reduction step.
    """
    has_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal has_nitro_reduction

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactant has nitro group and product has amine
                nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                try:
                    # For simplicity, assume single reactant for nitro reduction
                    if len(reactants) == 1:
                        reactant_mol = Chem.MolFromSmiles(reactants[0])
                        product_mol = Chem.MolFromSmiles(product)

                        if (
                            reactant_mol
                            and product_mol
                            and reactant_mol.HasSubstructMatch(nitro_pattern)
                            and product_mol.HasSubstructMatch(amine_pattern)
                            and not reactant_mol.HasSubstructMatch(amine_pattern)
                        ):
                            has_nitro_reduction = True
                            print("Nitro reduction detected")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Nitro reduction strategy: {has_nitro_reduction}")
    return has_nitro_reduction

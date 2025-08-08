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
    This function detects if the synthetic route involves reduction of a nitro group to an amine.
    """
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    amine_pattern = Chem.MolFromSmarts("[NH2]")

    nitro_to_amine = False

    def dfs_traverse(node):
        nonlocal nitro_to_amine

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Check for nitro in reactants
                    reactants_have_nitro = any(
                        Chem.MolFromSmiles(r)
                        and Chem.MolFromSmiles(r).HasSubstructMatch(nitro_pattern)
                        for r in reactants_smiles
                        if Chem.MolFromSmiles(r)
                    )

                    # Check for amine in product
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if (
                        reactants_have_nitro
                        and product_mol
                        and product_mol.HasSubstructMatch(amine_pattern)
                    ):
                        nitro_to_amine = True
                        print("Detected nitro to amine reduction")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitro_to_amine

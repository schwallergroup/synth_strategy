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
    This function detects a strategy involving benzylic bromination using NBS or similar reagents.
    """
    # Initialize tracking variable
    has_benzylic_bromination = False

    def dfs_traverse(node):
        nonlocal has_benzylic_bromination

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for benzylic bromination
            nbs_pattern = Chem.MolFromSmarts("O=C1CCC(=O)N1Br")
            benzylic_methyl_pattern = Chem.MolFromSmarts("c[CH3]")
            benzylic_bromomethyl_pattern = Chem.MolFromSmarts("c[CH2]Br")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            # Check if reactants include benzylic methyl and NBS
            has_benzylic_methyl = any(
                mol is not None and mol.HasSubstructMatch(benzylic_methyl_pattern)
                for mol in reactant_mols
            )
            has_nbs = any(
                mol is not None and mol.HasSubstructMatch(nbs_pattern) for mol in reactant_mols
            )

            # Check if product has benzylic bromomethyl
            product_has_bromomethyl = product_mol is not None and product_mol.HasSubstructMatch(
                benzylic_bromomethyl_pattern
            )

            if has_benzylic_methyl and has_nbs and product_has_bromomethyl:
                has_benzylic_bromination = True
                print("Detected benzylic bromination reaction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_benzylic_bromination

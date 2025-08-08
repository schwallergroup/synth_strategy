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
    This function detects pyrimidine ring formation in the synthetic route.
    """
    pyrimidine_formed = False

    def dfs_traverse(node):
        nonlocal pyrimidine_formed

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check if reactants don't have pyrimidine but product does
            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(r for r in reactants_mols):
                # SMARTS for pyrimidine
                pyrimidine_pattern = Chem.MolFromSmarts("c1ncncc1")

                # Check if product has pyrimidine
                product_has_pyrimidine = product_mol.HasSubstructMatch(pyrimidine_pattern)

                # Check if any reactant has pyrimidine
                reactants_have_pyrimidine = any(
                    r.HasSubstructMatch(pyrimidine_pattern) for r in reactants_mols if r
                )

                if product_has_pyrimidine and not reactants_have_pyrimidine:
                    print("Detected pyrimidine ring formation")
                    pyrimidine_formed = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return pyrimidine_formed

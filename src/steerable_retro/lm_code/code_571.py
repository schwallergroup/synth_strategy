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
    This function detects if the synthetic route involves pyrazole ring formation.
    """
    pyrazole_formed = False

    def dfs_traverse(node):
        nonlocal pyrazole_formed

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if reactants contain hydrazine derivative and product contains pyrazole
                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if product_mol:
                    # SMARTS for pyrazole
                    pyrazole_pattern = Chem.MolFromSmarts("[n]1[n][c][c][c]1")

                    # Check if product contains pyrazole
                    if product_mol.HasSubstructMatch(pyrazole_pattern):
                        # Check if any reactant has hydrazine pattern
                        hydrazine_pattern = Chem.MolFromSmarts("[NH][NH2]")
                        for reactant in reactants_mols:
                            if reactant and reactant.HasSubstructMatch(hydrazine_pattern):
                                print("Detected pyrazole formation from hydrazine derivative")
                                pyrazole_formed = True
                                break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return pyrazole_formed

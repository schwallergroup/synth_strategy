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
    This function detects multiple ester hydrolysis steps in a synthesis route.
    """
    ester_hydrolysis_count = 0

    def dfs_traverse(node):
        nonlocal ester_hydrolysis_count

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for ester pattern in reactants
            ester_pattern = Chem.MolFromSmarts("C(=O)OC")
            acid_pattern = Chem.MolFromSmarts("C(=O)O")
            methanol_pattern = Chem.MolFromSmarts("CO")

            has_ester = False
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(ester_pattern):
                        has_ester = True
                        break
                except:
                    continue

            # Check for acid pattern in product and methanol in reactants
            has_acid = False
            has_methanol = False

            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(acid_pattern):
                    has_acid = True
            except:
                pass

            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.GetNumAtoms() <= 3 and mol.HasSubstructMatch(methanol_pattern):
                        has_methanol = True
                except:
                    continue

            if has_ester and has_acid:
                print("Detected ester hydrolysis")
                ester_hydrolysis_count += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return ester_hydrolysis_count >= 2

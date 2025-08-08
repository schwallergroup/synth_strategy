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
    This function detects if the synthesis involves an N-oxide reduction step.
    """
    n_oxide_reduction = False

    def dfs_traverse(node):
        nonlocal n_oxide_reduction

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for N-oxide pattern in reactants
            n_oxide_pattern = Chem.MolFromSmarts("[n+][O-]")

            n_oxide_present = False
            for r_smi in reactants_smiles:
                try:
                    r_mol = Chem.MolFromSmiles(r_smi)
                    if r_mol and r_mol.HasSubstructMatch(n_oxide_pattern):
                        n_oxide_present = True
                        break
                except:
                    continue

            # Check if product has regular nitrogen where N-oxide was
            if n_oxide_present:
                try:
                    p_mol = Chem.MolFromSmiles(product_smiles)
                    if p_mol:
                        # If product doesn't have N-oxide but has pyridine, it's a reduction
                        if not p_mol.HasSubstructMatch(n_oxide_pattern) and p_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[n]")
                        ):
                            n_oxide_reduction = True
                            print("Detected N-oxide reduction")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return n_oxide_reduction

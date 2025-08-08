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
    This function detects if the synthetic route involves multiple SNAr reactions
    for fragment assembly.
    """
    snar_reactions_count = 0

    def dfs_traverse(node):
        nonlocal snar_reactions_count

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for SNAr pattern: halogen on electron-deficient aromatic + nucleophile
            # Look for patterns like F/Cl on aromatic ring being replaced by N nucleophile
            try:
                product_mol = Chem.MolFromSmiles(product)

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    # Check for halogen on aromatic
                    halo_aromatic_pattern = Chem.MolFromSmarts(
                        "[F,Cl,Br,I]c1[c,n][c,n][c,n][c,n][c,n]1"
                    )

                    # Check for nitrogen nucleophile
                    n_nucleophile_pattern = Chem.MolFromSmarts("[NH,NH2][c,C]")

                    if reactant_mol.HasSubstructMatch(halo_aromatic_pattern) and any(
                        Chem.MolFromSmiles(r).HasSubstructMatch(n_nucleophile_pattern)
                        for r in reactants
                        if Chem.MolFromSmiles(r)
                    ):

                        # Check if product has C-N bond where halogen was
                        if product_mol and product_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("Nc1[c,n][c,n][c,n][c,n][c,n]1")
                        ):
                            snar_reactions_count += 1
                            print(f"Detected SNAr reaction: {rsmi}")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return snar_reactions_count >= 2

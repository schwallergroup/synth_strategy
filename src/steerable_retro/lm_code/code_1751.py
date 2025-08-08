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
    This function detects formation of heterocyclic rings, particularly imidazopyridine.
    """
    heterocycle_formation_found = False

    def dfs_traverse(node):
        nonlocal heterocycle_formation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Count rings in reactants and product
            try:
                p_mol = Chem.MolFromSmiles(product)
                if p_mol:
                    product_rings = Chem.GetSSSR(p_mol)

                    total_reactant_rings = 0
                    for reactant in reactants:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol:
                            total_reactant_rings += len(Chem.GetSSSR(r_mol))

                    # Check if product has more rings than reactants combined
                    if len(product_rings) > total_reactant_rings:
                        # Check specifically for imidazopyridine pattern
                        imidazopyridine_pattern = Chem.MolFromSmarts(
                            "[n]1[c][n][c]2[c]1[c][c][c][n]2"
                        )
                        if p_mol.HasSubstructMatch(imidazopyridine_pattern):
                            print("Found imidazopyridine formation")
                            heterocycle_formation_found = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return heterocycle_formation_found

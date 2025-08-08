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
    This function detects a synthetic strategy involving thioether formation
    between an aromatic thiol and an alkyl halide.
    """
    thioether_found = False

    def dfs_traverse(node):
        nonlocal thioether_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a thioether formation reaction
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Look for thioether pattern in product
                    thioether_patt = Chem.MolFromSmarts("[#6]-[#16]-[#6]")
                    if product_mol.HasSubstructMatch(thioether_patt):
                        # Check if reactants contain aromatic thiol and alkyl halide
                        aromatic_thiol_patt = Chem.MolFromSmarts("[c]-[#16H]")
                        alkyl_halide_patt = Chem.MolFromSmarts("[#6]-[#17,#35,#53]")

                        has_aromatic_thiol = False
                        has_alkyl_halide = False

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                if reactant_mol.HasSubstructMatch(aromatic_thiol_patt):
                                    has_aromatic_thiol = True
                                if reactant_mol.HasSubstructMatch(alkyl_halide_patt):
                                    has_alkyl_halide = True

                        if has_aromatic_thiol and has_alkyl_halide:
                            print("Found thioether formation reaction")
                            thioether_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return thioether_found

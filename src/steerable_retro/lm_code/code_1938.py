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
    This function detects borylation reactions that prepare coupling partners.
    It looks for reactions converting aryl halides to boronic esters.
    """
    borylation_found = False

    def dfs_traverse(node):
        nonlocal borylation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl halide pattern in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")

            # Check for boronic ester pattern in product
            boronic_ester_pattern = Chem.MolFromSmarts("[c][B]1[O][C]([C])([C])[C]([C])([C])[O]1")

            has_aryl_halide = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True
                        break
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product)
                has_boronic_ester = product_mol and product_mol.HasSubstructMatch(
                    boronic_ester_pattern
                )
            except:
                has_boronic_ester = False

            if has_aryl_halide and has_boronic_ester:
                print(f"Found borylation reaction: {rsmi}")
                borylation_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Borylation preparation found: {borylation_found}")
    return borylation_found

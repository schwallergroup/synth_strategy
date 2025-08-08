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
    Detects if the synthesis involves N-sulfonylation of a pyrazole ring.
    """
    found_sulfonylation = False

    def dfs_traverse(node):
        nonlocal found_sulfonylation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            products = rsmi.split(">")[-1]

            try:
                # Pattern for pyrazole NH
                pyrazole_nh_pattern = Chem.MolFromSmarts("[nH]1cncc1")
                # Pattern for sulfonyl chloride
                sulfonyl_cl_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])-[Cl]")
                # Pattern for sulfonylated pyrazole
                sulfonylated_pyrazole_pattern = Chem.MolFromSmarts(
                    "[#16](=[#8])(=[#8])-[n]1[c][n][c][c]1"
                )

                pyrazole_found = False
                sulfonyl_cl_found = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(pyrazole_nh_pattern):
                            pyrazole_found = True
                        if reactant_mol.HasSubstructMatch(sulfonyl_cl_pattern):
                            sulfonyl_cl_found = True

                product_mol = Chem.MolFromSmiles(products)
                sulfonylated_product = False
                if product_mol and product_mol.HasSubstructMatch(sulfonylated_pyrazole_pattern):
                    sulfonylated_product = True

                if pyrazole_found and sulfonyl_cl_found and sulfonylated_product:
                    found_sulfonylation = True
                    print("Found pyrazole N-sulfonylation reaction")
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_sulfonylation

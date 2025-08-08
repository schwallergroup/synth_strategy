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
    This function detects the formation of a benzimidazole ring system through
    the coupling of a diamine with an isothiocyanate.
    """
    found_diamine_isothiocyanate_coupling = False

    def dfs_traverse(node):
        nonlocal found_diamine_isothiocyanate_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains benzimidazole
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                benzimidazole_pattern = Chem.MolFromSmarts("[nH]1cnc2ccccc12")
                if product_mol.HasSubstructMatch(benzimidazole_pattern):
                    # Check if reactants contain diamine and isothiocyanate
                    diamine_found = False
                    isothiocyanate_found = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            # Check for diamine (two NH2 groups)
                            diamine_pattern = Chem.MolFromSmarts("[NH2].[NH2]")
                            if reactant_mol.HasSubstructMatch(diamine_pattern):
                                diamine_found = True

                            # Check for isothiocyanate (N=C=S)
                            isothiocyanate_pattern = Chem.MolFromSmarts("[#6][N]=[C]=[S]")
                            if reactant_mol.HasSubstructMatch(isothiocyanate_pattern):
                                isothiocyanate_found = True

                    if diamine_found and isothiocyanate_found:
                        found_diamine_isothiocyanate_coupling = True
                        print("Found benzimidazole formation via diamine-isothiocyanate coupling")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_diamine_isothiocyanate_coupling

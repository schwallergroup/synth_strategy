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
    This function detects if the synthetic route involves coupling a pyridine with a dichlorophenyl group.
    """
    pyridine_dichlorophenyl_coupling_found = False

    def dfs_traverse(node):
        nonlocal pyridine_dichlorophenyl_coupling_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check reactants for pyridine and dichlorophenyl
            pyridine_found = False
            dichlorophenyl_found = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if not reactant_mol:
                    continue

                # Check for pyridine
                pyridine_pattern = Chem.MolFromSmarts("[n]1[c][c][c][c][c]1")
                if reactant_mol.HasSubstructMatch(pyridine_pattern):
                    pyridine_found = True

                # Check for dichlorophenyl
                dichlorophenyl_pattern = Chem.MolFromSmarts("[c]1([Cl])[c]([Cl])[c][c][c][c]1")
                if reactant_mol.HasSubstructMatch(dichlorophenyl_pattern):
                    dichlorophenyl_found = True

            # Check product for coupled structure
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                pyridine_pattern = Chem.MolFromSmarts("[n]1[c][c][c][c][c]1")
                dichlorophenyl_pattern = Chem.MolFromSmarts("[c]1([Cl])[c]([Cl])[c][c][c][c]1")

                if (
                    pyridine_found
                    and dichlorophenyl_found
                    and product_mol.HasSubstructMatch(pyridine_pattern)
                    and product_mol.HasSubstructMatch(dichlorophenyl_pattern)
                ):
                    print("Detected pyridine-dichlorophenyl coupling")
                    pyridine_dichlorophenyl_coupling_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return pyridine_dichlorophenyl_coupling_found

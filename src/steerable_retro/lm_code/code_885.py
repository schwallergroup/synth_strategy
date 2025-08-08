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
    This function detects Suzuki coupling for biaryl formation.
    Looks for C-C bond formation between two aromatic rings where one reactant
    contains a boronic acid and the other is an activated arene.
    """
    suzuki_detected = False

    def dfs_traverse(node):
        nonlocal suzuki_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronic acid in reactants
                boronic_acid_pattern = Chem.MolFromSmarts("[c]-[B]([OH])[OH]")

                # Check for triflate or halide in reactants
                leaving_group_pattern = Chem.MolFromSmarts("[c]-[O,F,Cl,Br,I]")

                # Check for biaryl in product
                biaryl_pattern = Chem.MolFromSmarts("c1ccccc1-c1ccccc1")

                has_boronic_acid = False
                has_leaving_group = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(boronic_acid_pattern):
                            has_boronic_acid = True
                        if mol and mol.HasSubstructMatch(leaving_group_pattern):
                            has_leaving_group = True
                    except:
                        continue

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    has_biaryl = product_mol and product_mol.HasSubstructMatch(biaryl_pattern)

                    if has_boronic_acid and has_leaving_group and has_biaryl:
                        print("Detected Suzuki coupling for biaryl formation")
                        suzuki_detected = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_detected

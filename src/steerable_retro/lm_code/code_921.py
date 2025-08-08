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
    Detects if the synthesis includes a Suzuki coupling with a cyclopropylboronic acid.
    """
    found = False

    def dfs_traverse(node):
        nonlocal found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl bromide pattern
            aryl_bromide_pattern = Chem.MolFromSmarts("c[Br]")
            # Check for cyclopropylboronic acid or derivative
            cyclopropyl_pattern = Chem.MolFromSmarts("[C]1[C][C]1")
            boron_pattern = Chem.MolFromSmarts("[B]")

            has_aryl_bromide = False
            has_cyclopropyl_boron = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(aryl_bromide_pattern):
                        has_aryl_bromide = True
                    if mol.HasSubstructMatch(cyclopropyl_pattern) and mol.HasSubstructMatch(
                        boron_pattern
                    ):
                        has_cyclopropyl_boron = True

            if has_aryl_bromide and has_cyclopropyl_boron:
                # Check if product has cyclopropyl attached to aryl
                cyclopropyl_aryl_pattern = Chem.MolFromSmarts("c[C]1[C][C]1")
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol and prod_mol.HasSubstructMatch(cyclopropyl_aryl_pattern):
                    print("Found Suzuki coupling with cyclopropyl group")
                    found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found

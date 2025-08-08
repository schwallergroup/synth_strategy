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
    Detects if the synthesis includes a Suzuki coupling reaction.
    """
    has_suzuki_coupling = False

    def dfs_traverse(node):
        nonlocal has_suzuki_coupling

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Look for boronic acid in reactants
            has_boronic_acid = False
            has_aryl_halide = False

            for reactant in reactants:
                if "B(O)O" in reactant or "OB(O)" in reactant:
                    has_boronic_acid = True

                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    aryl_halide_pattern = Chem.MolFromSmarts("c[Cl,Br,I]")
                    if reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True

            # Check if product has biaryl system
            if has_boronic_acid and has_aryl_halide:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    biaryl_pattern = Chem.MolFromSmarts("c:c-c:c")
                    if product_mol.HasSubstructMatch(biaryl_pattern):
                        has_suzuki_coupling = True
                        print("Detected Suzuki coupling")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_suzuki_coupling

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
    This function detects Suzuki coupling reactions (aryl-halide + boronic acid/ester).
    """
    has_suzuki_coupling = False

    def dfs_traverse(node):
        nonlocal has_suzuki_coupling

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aryl halide pattern
                aryl_halide_pattern = Chem.MolFromSmarts("c[Br,I,Cl]")

                # Check for boronic acid/ester pattern
                boronic_pattern = Chem.MolFromSmarts("cB(O[#6])(O[#6])")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                has_aryl_halide = any(
                    r and r.HasSubstructMatch(aryl_halide_pattern) for r in reactant_mols
                )
                has_boronic = any(r and r.HasSubstructMatch(boronic_pattern) for r in reactant_mols)

                # If both patterns are found in reactants and product has a new C-C bond
                if has_aryl_halide and has_boronic and product_mol:
                    has_suzuki_coupling = True
                    print("Detected Suzuki coupling reaction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    return has_suzuki_coupling

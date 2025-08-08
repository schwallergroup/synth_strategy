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
    Detects if the route contains Suzuki coupling reactions forming biaryl bonds.
    Looks for reactions where a boronic acid/ester and aryl halide form a biaryl C-C bond.
    """
    suzuki_count = 0

    def dfs_traverse(node):
        nonlocal suzuki_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for boronic acid/ester pattern in reactants
            boronic_pattern = Chem.MolFromSmarts("[c,C]-[B]([O,OH])[O,OH]")
            boronic_ester_pattern = Chem.MolFromSmarts("[c,C]-[B]1O[C][C]O1")

            # Check for aryl halide pattern in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Cl,Br,I]")

            # Check for biaryl pattern in product
            biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")

            has_boronic = False
            has_aryl_halide = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(boronic_pattern) or mol.HasSubstructMatch(
                            boronic_ester_pattern
                        ):
                            has_boronic = True
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product)
                if (
                    product_mol
                    and product_mol.HasSubstructMatch(biaryl_pattern)
                    and has_boronic
                    and has_aryl_halide
                ):
                    suzuki_count += 1
                    print(
                        f"Found Suzuki coupling forming biaryl bond at depth {node.get('depth', 'unknown')}"
                    )
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total Suzuki couplings forming biaryl bonds: {suzuki_count}")
    return suzuki_count >= 1

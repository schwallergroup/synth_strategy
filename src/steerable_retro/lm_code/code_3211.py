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
    Detects a strategy involving construction of multiple heterocycles
    with at least one heterocycle formed in the first half of synthesis
    and another heterocycle formed in the second half.
    """
    heterocycle_formations = []
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for heterocycle formation
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants.split(".")]
            product_mol = Chem.MolFromSmiles(product)

            if product_mol and all(r for r in reactant_mols):
                # Count rings in reactants and product
                reactant_rings = sum(len(Chem.GetSSSR(r)) for r in reactant_mols if r)
                product_rings = len(Chem.GetSSSR(product_mol)) if product_mol else 0

                # If product has more rings than reactants combined, a ring was formed
                if product_rings > reactant_rings:
                    # Check if the new ring is a heterocycle
                    imidazole_pattern = Chem.MolFromSmarts("[n]1[c][n][c][c]1")
                    thiazole_pattern = Chem.MolFromSmarts("[s]1[c][n][c][c]1")
                    pyrrole_pattern = Chem.MolFromSmarts("[nH]1[c][c][c][c]1")

                    # Check for various heterocycles in product but not in reactants
                    product_has_imidazole = product_mol.HasSubstructMatch(imidazole_pattern)
                    product_has_thiazole = product_mol.HasSubstructMatch(thiazole_pattern)
                    product_has_pyrrole = product_mol.HasSubstructMatch(pyrrole_pattern)

                    reactants_have_imidazole = any(
                        r and r.HasSubstructMatch(imidazole_pattern) for r in reactant_mols
                    )
                    reactants_have_thiazole = any(
                        r and r.HasSubstructMatch(thiazole_pattern) for r in reactant_mols
                    )
                    reactants_have_pyrrole = any(
                        r and r.HasSubstructMatch(pyrrole_pattern) for r in reactant_mols
                    )

                    if (
                        (product_has_imidazole and not reactants_have_imidazole)
                        or (product_has_thiazole and not reactants_have_thiazole)
                        or (product_has_pyrrole and not reactants_have_pyrrole)
                    ):
                        heterocycle_formations.append((depth, "heterocycle_formation"))
                        print(f"Heterocycle formation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have heterocycle formations in both early and late stages
    if not heterocycle_formations:
        return False

    # Sort by depth
    heterocycle_formations.sort(key=lambda x: x[0])

    # Calculate midpoint of synthesis
    midpoint = max_depth / 2

    early_formations = [f for f in heterocycle_formations if f[0] > midpoint]
    late_formations = [f for f in heterocycle_formations if f[0] <= midpoint]

    # Check if we have heterocycle formations in both early and late stages
    has_multi_heterocycle_strategy = len(early_formations) > 0 and len(late_formations) > 0

    print(f"Multi-heterocycle construction strategy detected: {has_multi_heterocycle_strategy}")
    print(f"Early formations: {early_formations}")
    print(f"Late formations: {late_formations}")

    return has_multi_heterocycle_strategy

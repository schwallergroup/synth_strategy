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
    This function detects if the synthetic route includes early aromatic functionalization
    (nitration, halogenation) in the first half of the synthesis.
    """
    functionalization_depths = []
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal functionalization_depths, max_depth

        # Track maximum depth to determine early vs late
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Patterns for nitration and halogenation
            nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")
            bromo_pattern = Chem.MolFromSmarts("c[Br]")
            fluoro_pattern = Chem.MolFromSmarts("c[F]")

            # Check if product has new functional group that reactants don't
            product_mol = None
            try:
                product_mol = Chem.MolFromSmiles(product)
            except:
                pass

            if product_mol:
                # Check for nitration
                if product_mol.HasSubstructMatch(nitro_pattern):
                    reactant_has_nitro = False
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(nitro_pattern):
                                reactant_has_nitro = True
                                break
                        except:
                            continue

                    if not reactant_has_nitro:
                        print(f"Detected nitration at depth {depth}")
                        functionalization_depths.append(depth)

                # Check for bromination
                if product_mol.HasSubstructMatch(bromo_pattern):
                    reactant_has_bromo = False
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(bromo_pattern):
                                reactant_has_bromo = True
                                break
                        except:
                            continue

                    if not reactant_has_bromo:
                        print(f"Detected bromination at depth {depth}")
                        functionalization_depths.append(depth)

                # Check for fluorination
                if product_mol.HasSubstructMatch(fluoro_pattern):
                    reactant_has_fluoro = False
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(fluoro_pattern):
                                reactant_has_fluoro = True
                                break
                        except:
                            continue

                    if not reactant_has_fluoro:
                        print(f"Detected fluorination at depth {depth}")
                        functionalization_depths.append(depth)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if functionalization occurs in the second half of the synthesis (higher depth)
    # Since we're traversing retrosynthetically, higher depth means earlier in the actual synthesis
    if functionalization_depths and max_depth > 0:
        # Consider it early if the average depth is in the second half
        avg_depth = sum(functionalization_depths) / len(functionalization_depths)
        return avg_depth > (max_depth / 2)

    return False

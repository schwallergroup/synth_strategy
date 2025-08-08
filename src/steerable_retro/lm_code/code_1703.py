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
    Detects a synthetic strategy where a piperazine is introduced via SNAr in the late stage of synthesis.
    Looks for:
    1. SNAr reaction (halogen displacement by nitrogen nucleophile)
    2. Piperazine as the nucleophile
    3. Reaction occurring in the late stage (low depth)
    4. Activated aromatic ring (with electron-withdrawing groups)
    """
    # Track if we found the SNAr with piperazine
    found_snar_with_piperazine = False

    def dfs_traverse(node, depth=0):
        nonlocal found_snar_with_piperazine

        if node["type"] == "reaction" and depth <= 2:  # Focus on late-stage reactions (low depth)
            # Extract reactants and product
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is potentially an SNAr reaction
                product_mol = Chem.MolFromSmiles(product_smiles)

                # Look for piperazine in reactants
                piperazine_pattern = Chem.MolFromSmarts(
                    "[N]1CCN([C])CC1"
                )  # Methylpiperazine pattern

                has_piperazine = False
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(piperazine_pattern):
                        has_piperazine = True
                        break

                # Look for halogen in other reactant (potential leaving group)
                halogen_pattern = Chem.MolFromSmarts("c[F,Cl,Br,I]")

                has_aryl_halide = False
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(halogen_pattern):
                        has_aryl_halide = True
                        break

                # Check if product has piperazine attached to aromatic ring
                piperazine_on_aryl_pattern = Chem.MolFromSmarts("c[N]1CCN([C])CC1")

                if (
                    has_piperazine
                    and has_aryl_halide
                    and product_mol
                    and product_mol.HasSubstructMatch(piperazine_on_aryl_pattern)
                ):
                    print(f"Found SNAr with piperazine at depth {depth}")

                    # Check for activating groups (nitro, fluoro) on the aromatic ring
                    nitro_pattern = Chem.MolFromSmarts("c[N+](=O)[O-]")
                    fluoro_pattern = Chem.MolFromSmarts("cF")

                    has_activating_groups = False
                    for reactant in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(
                                nitro_pattern
                            ) or reactant_mol.HasSubstructMatch(fluoro_pattern):
                                has_activating_groups = True
                                break

                    if has_activating_groups:
                        print("Found activating groups for SNAr")
                        found_snar_with_piperazine = True

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_snar_with_piperazine

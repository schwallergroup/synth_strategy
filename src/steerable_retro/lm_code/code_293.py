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
    This function detects a sequence where a nitro group is reduced to an amine
    and later converted to a sulfonamide.
    """
    # Track reactions with their depths
    reactions_by_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                reactions_by_depth[depth] = node["metadata"]["rsmi"]

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal to collect reactions
    dfs_traverse(route)

    # Check for nitro reduction
    found_nitro_reduction = False
    nitro_reduction_depth = None

    # Check for sulfonamide formation
    found_sulfonamide_formation = False
    sulfonamide_formation_depth = None

    # Analyze reactions
    for depth, rsmi in reactions_by_depth.items():
        reactants = rsmi.split(">")[0].split(".")
        product = rsmi.split(">")[-1]

        # Check for nitro reduction
        product_mol = Chem.MolFromSmiles(product) if product else None
        if product_mol:
            amine_pattern = Chem.MolFromSmarts("[#7;H2]")
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")

            reactants_have_nitro = False
            for mol in reactant_mols:
                if mol and mol.HasSubstructMatch(nitro_pattern):
                    reactants_have_nitro = True
                    break

            if product_mol.HasSubstructMatch(amine_pattern) and reactants_have_nitro:
                found_nitro_reduction = True
                nitro_reduction_depth = depth
                print(f"Found nitro reduction at depth {depth}")

        # Check for sulfonamide formation
        product_mol = Chem.MolFromSmiles(product) if product else None
        if product_mol:
            sulfonamide_pattern = Chem.MolFromSmarts("[#7]S(=O)(=O)[#6]")
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            amine_pattern = Chem.MolFromSmarts("[#7;H2]")

            reactants_have_amine = False
            for mol in reactant_mols:
                if mol and mol.HasSubstructMatch(amine_pattern):
                    reactants_have_amine = True
                    break

            if product_mol.HasSubstructMatch(sulfonamide_pattern) and reactants_have_amine:
                found_sulfonamide_formation = True
                sulfonamide_formation_depth = depth
                print(f"Found sulfonamide formation at depth {depth}")

    # Check if we found both transformations in the correct order
    if (
        found_nitro_reduction
        and found_sulfonamide_formation
        and nitro_reduction_depth is not None
        and sulfonamide_formation_depth is not None
        and sulfonamide_formation_depth < nitro_reduction_depth
    ):
        print(
            f"Found nitro-to-sulfonamide sequence: nitro reduction at depth {nitro_reduction_depth}, sulfonamide formation at depth {sulfonamide_formation_depth}"
        )
        return True

    return False

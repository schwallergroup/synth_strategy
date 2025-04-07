#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects a strategy where an olefination reaction is used
    to form a C=C bond, followed by an amide bond formation.
    """
    olefination_depth = None
    amide_formation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal olefination_depth, amide_formation_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for olefination reaction (phosphonate ester + carbonyl)
                phosphonate_pattern = Chem.MolFromSmarts("[P](=[O])([O][C])[O][C]")
                carbonyl_pattern = Chem.MolFromSmarts("[C](=[O])[C,H]")
                alkene_pattern = Chem.MolFromSmarts("[C]=[C]")

                has_phosphonate = False
                has_carbonyl = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(phosphonate_pattern):
                            has_phosphonate = True
                        if mol and mol.HasSubstructMatch(carbonyl_pattern):
                            has_carbonyl = True
                    except:
                        continue

                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if (
                        has_phosphonate
                        and has_carbonyl
                        and prod_mol
                        and prod_mol.HasSubstructMatch(alkene_pattern)
                    ):
                        print(f"Found olefination reaction at depth {depth}")
                        olefination_depth = depth
                except:
                    pass

                # Check for amide formation
                carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                amide_pattern = Chem.MolFromSmarts("[C](=[O])[NH]")

                has_acid = False
                has_amine = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(carboxylic_acid_pattern):
                            has_acid = True
                        if mol and mol.HasSubstructMatch(amine_pattern):
                            has_amine = True
                    except:
                        continue

                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if (
                        has_acid
                        and has_amine
                        and prod_mol
                        and prod_mol.HasSubstructMatch(amide_pattern)
                    ):
                        print(f"Found amide formation at depth {depth}")
                        amide_formation_depth = depth
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if olefination occurs before amide formation
    if olefination_depth is not None and amide_formation_depth is not None:
        if (
            olefination_depth < amide_formation_depth
        ):  # Remember: lower depth = later in synthesis
            print(
                f"Confirmed olefination-amide strategy: olefination at {olefination_depth}, amide formation at {amide_formation_depth}"
            )
            return True
    return False

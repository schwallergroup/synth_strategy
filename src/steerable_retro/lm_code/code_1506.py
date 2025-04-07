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
    This function detects a linear synthesis strategy with sequential
    functional group transformations (ester → aldehyde → imine → hydroxamic amide).
    """
    # Track functional group transformations
    transformations = []

    def dfs_traverse(node):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check number of reactants (linear synthesis typically has 1-2 reactants)
            if len(reactants) <= 2:
                # Identify functional groups
                r_mol = Chem.MolFromSmiles(reactants[0]) if reactants else None
                p_mol = Chem.MolFromSmiles(product) if product else None

                if r_mol and p_mol:
                    # Check for ester
                    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CH3]")
                    # Check for aldehyde
                    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=[OX1])")
                    # Check for imine/oxime
                    imine_pattern = Chem.MolFromSmarts("[NX2]=[CX3]")
                    # Check for hydroxamic amide
                    hydroxamic_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[NX3][OX2]")

                    # Record transformations
                    if r_mol.HasSubstructMatch(ester_pattern) and not p_mol.HasSubstructMatch(
                        ester_pattern
                    ):
                        transformations.append("ester_removal")
                        print("Found ester removal transformation")

                    if not r_mol.HasSubstructMatch(aldehyde_pattern) and p_mol.HasSubstructMatch(
                        aldehyde_pattern
                    ):
                        transformations.append("aldehyde_formation")
                        print("Found aldehyde formation")

                    if r_mol.HasSubstructMatch(aldehyde_pattern) and not p_mol.HasSubstructMatch(
                        aldehyde_pattern
                    ):
                        transformations.append("aldehyde_transformation")
                        print("Found aldehyde transformation")

                    if not r_mol.HasSubstructMatch(imine_pattern) and p_mol.HasSubstructMatch(
                        imine_pattern
                    ):
                        transformations.append("imine_formation")
                        print("Found imine formation")

                    if r_mol.HasSubstructMatch(imine_pattern) and not p_mol.HasSubstructMatch(
                        imine_pattern
                    ):
                        transformations.append("imine_transformation")
                        print("Found imine transformation")

                    if not r_mol.HasSubstructMatch(hydroxamic_pattern) and p_mol.HasSubstructMatch(
                        hydroxamic_pattern
                    ):
                        transformations.append("hydroxamic_formation")
                        print("Found hydroxamic amide formation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have the expected sequence of transformations
    # We need at least 3 transformations to consider this a significant FG transformation strategy
    has_expected_sequence = len(transformations) >= 3

    # Check if we have both formation and transformation of key intermediates
    has_aldehyde_cycle = (
        "aldehyde_formation" in transformations and "aldehyde_transformation" in transformations
    )
    has_imine_cycle = (
        "imine_formation" in transformations and "imine_transformation" in transformations
    )
    has_hydroxamic = "hydroxamic_formation" in transformations

    return has_expected_sequence and (has_aldehyde_cycle or has_imine_cycle) and has_hydroxamic

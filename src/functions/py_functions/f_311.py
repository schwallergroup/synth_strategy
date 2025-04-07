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
    This function detects aromatic C-N bond formation via SNAr or coupling reactions.
    """
    aromatic_c_n_formation_detected = False

    def dfs_traverse(node):
        nonlocal aromatic_c_n_formation_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            try:
                reactants = reactants_part.split(".")
                product_mol = Chem.MolFromSmiles(product_part)

                # Aromatic C-N bond pattern
                aromatic_c_n_pattern = Chem.MolFromSmarts("[c][N]")

                # Aromatic halide pattern (for SNAr or coupling precursor)
                aromatic_halide_pattern = Chem.MolFromSmarts("[c][Br,Cl,I,F]")

                # Amine pattern
                amine_pattern = Chem.MolFromSmarts("[N;H1,H2]")

                # Check if product has aromatic C-N bond
                if product_mol and product_mol.HasSubstructMatch(aromatic_c_n_pattern):
                    # Check if reactants have aromatic halide and amine
                    has_aromatic_halide = False
                    has_amine = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if not reactant_mol:
                            continue

                        if reactant_mol.HasSubstructMatch(aromatic_halide_pattern):
                            has_aromatic_halide = True

                        if reactant_mol.HasSubstructMatch(amine_pattern):
                            has_amine = True

                    if has_aromatic_halide and has_amine:
                        aromatic_c_n_formation_detected = True
                        print("Detected aromatic C-N bond formation")
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return aromatic_c_n_formation_detected

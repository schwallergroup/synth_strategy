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
    This function detects thioether formation via nucleophilic substitution.
    """
    thioether_formation_found = False

    def dfs_traverse(node):
        nonlocal thioether_formation_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for thioether formation
                thiol_pattern = Chem.MolFromSmarts("[SH]C")
                alkyl_halide_pattern = Chem.MolFromSmarts("*CCC[Br,Cl,I,F]")
                thioether_pattern = Chem.MolFromSmarts("*CCCS(C)")

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if (
                        product_mol
                        and thiol_pattern
                        and alkyl_halide_pattern
                        and thioether_pattern
                    ):
                        thiol_found = False
                        alkyl_halide_found = False

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                if reactant_mol.HasSubstructMatch(thiol_pattern):
                                    thiol_found = True
                                if reactant_mol.HasSubstructMatch(alkyl_halide_pattern):
                                    alkyl_halide_found = True

                        if (
                            thiol_found
                            and alkyl_halide_found
                            and product_mol.HasSubstructMatch(thioether_pattern)
                        ):
                            print(
                                "Found thioether formation via nucleophilic substitution"
                            )
                            thioether_formation_found = True
                except:
                    print(
                        "Error in SMILES processing for thioether formation detection"
                    )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return thioether_formation_found

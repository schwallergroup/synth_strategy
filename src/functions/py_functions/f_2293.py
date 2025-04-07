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
    This function detects late-stage thiazole formation via S-alkylation.
    """
    thiazole_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal thiazole_formation_detected

        if node["type"] == "reaction" and depth <= 1:  # Late stage (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains thiazole
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("c1sccn1")
                ):
                    # Check if reactants contain thiol and alkyl halide
                    for reactant in reactants_smiles:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(
                                Chem.MolFromSmarts("[SH]")
                            ):
                                for r2 in reactants_smiles:
                                    r2_mol = Chem.MolFromSmiles(r2)
                                    if r2_mol and r2_mol.HasSubstructMatch(
                                        Chem.MolFromSmarts("[C][Cl,Br,I]")
                                    ):
                                        thiazole_formation_detected = True
                                        print(
                                            f"Late-stage thiazole formation detected at depth {depth}"
                                        )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return thiazole_formation_detected

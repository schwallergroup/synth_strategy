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
    This function detects pyrazole ring formation from hydrazine and carbonyl compounds.
    """
    pyrazole_formation_detected = False

    def dfs_traverse(node):
        nonlocal pyrazole_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if hydrazine is in reactants
                hydrazine_pattern = Chem.MolFromSmarts("[NH2][NH2]")
                carbonyl_pattern = Chem.MolFromSmarts("[#6][#6](=[O])[#6,#1]")
                pyrazole_pattern = Chem.MolFromSmarts("[#7]1[#7][#6][#6][#6]1")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol and any(
                    r and r.HasSubstructMatch(hydrazine_pattern) for r in reactant_mols
                ):
                    if any(
                        r and r.HasSubstructMatch(carbonyl_pattern)
                        for r in reactant_mols
                    ):
                        if product_mol.HasSubstructMatch(pyrazole_pattern):
                            print(
                                "Detected pyrazole formation from hydrazine and carbonyl compounds"
                            )
                            pyrazole_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return pyrazole_formation_detected

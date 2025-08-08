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
    Detects a strategy involving the use of a perfluorinated phenol as a leaving group
    in sulfonamide formation.
    """
    # Track if we find the perfluorinated leaving group strategy
    found_perfluorophenol = False
    found_perfluorophenyl_sulfonate = False
    found_sulfonamide_formation = False

    # SMARTS patterns
    perfluorophenol_pattern = Chem.MolFromSmarts("[OH][c]1[c]([F])[c]([F])[c]([F])[c]([F])[c]1[F]")
    perfluorophenyl_sulfonate_pattern = Chem.MolFromSmarts(
        "[#16](=[#8])(=[#8])[#8][c]1[c]([F])[c]([F])[c]([F])[c]([F])[c]1[F]"
    )

    def dfs_traverse(node):
        nonlocal found_perfluorophenol, found_perfluorophenyl_sulfonate, found_sulfonamide_formation

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if mol.HasSubstructMatch(perfluorophenol_pattern):
                    found_perfluorophenol = True
                if mol.HasSubstructMatch(perfluorophenyl_sulfonate_pattern):
                    found_perfluorophenyl_sulfonate = True

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this reaction converts a perfluorophenyl sulfonate to a sulfonamide
            has_perfluorophenyl_sulfonate = False
            has_amine = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(perfluorophenyl_sulfonate_pattern):
                        has_perfluorophenyl_sulfonate = True
                    if mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")):
                        has_amine = True

            product_mol = Chem.MolFromSmiles(product)
            has_sulfonamide = product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts("[#16](=[#8])(=[#8])[#7]")
            )

            if has_perfluorophenyl_sulfonate and has_amine and has_sulfonamide:
                found_sulfonamide_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we found the complete strategy
    result = (
        found_perfluorophenol and found_perfluorophenyl_sulfonate and found_sulfonamide_formation
    )

    print(f"Perfluorinated leaving group strategy detected: {result}")
    print(f"  - Found perfluorophenol: {found_perfluorophenol}")
    print(f"  - Found perfluorophenyl sulfonate: {found_perfluorophenyl_sulfonate}")
    print(f"  - Found sulfonamide formation: {found_sulfonamide_formation}")

    return result

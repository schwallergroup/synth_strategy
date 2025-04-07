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
    This function detects a synthetic strategy involving sequential functionalization
    of pyridine rings, including N-oxide formation and aromatic substitution.
    """
    # Track if we found each key step
    pyridine_present = False
    n_oxide_formation = False
    aromatic_substitution = False

    def dfs_traverse(node):
        nonlocal pyridine_present, n_oxide_formation, aromatic_substitution

        if node["type"] == "mol":
            # Check for pyridine in any molecule
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                pyridine_pattern = Chem.MolFromSmarts("c1ccccn1")
                if mol.HasSubstructMatch(pyridine_pattern):
                    print("Found pyridine structure")
                    pyridine_present = True

        elif node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check for N-oxide formation
            if pyridine_present and not n_oxide_formation:
                product_mol = Chem.MolFromSmiles(product_part)
                n_oxide_pattern = Chem.MolFromSmarts("[n+]-[O-]")
                if product_mol and product_mol.HasSubstructMatch(n_oxide_pattern):
                    print("Found N-oxide formation")
                    n_oxide_formation = True

            # Check for aromatic substitution (any new group on aromatic ring)
            if n_oxide_formation and not aromatic_substitution:
                reactants = reactants_part.split(".")
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product_part)

                # Look for nitration, halogenation, or other aromatic substitutions
                nitro_pattern = Chem.MolFromSmarts("[#6]~[N+](~[O-])~[O-]")
                halogen_pattern = Chem.MolFromSmarts("[#6]-[F,Cl,Br,I]")
                ether_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]")

                if product_mol:
                    if (
                        product_mol.HasSubstructMatch(nitro_pattern)
                        or product_mol.HasSubstructMatch(halogen_pattern)
                        or product_mol.HasSubstructMatch(ether_pattern)
                    ):
                        print("Found aromatic substitution")
                        aromatic_substitution = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if all steps of the strategy were found
    return pyridine_present and n_oxide_formation and aromatic_substitution

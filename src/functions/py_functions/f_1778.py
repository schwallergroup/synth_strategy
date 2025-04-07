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
    This function detects a synthetic strategy involving N-alkylation
    of a secondary amine to form a tertiary amine.
    """
    has_n_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_n_alkylation

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for secondary amine in reactants
                has_secondary_amine = False
                for reactant in reactants_smiles:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[N;H1]")):
                        has_secondary_amine = True
                        break

                # Check for benzyl or alkyl group in reactants
                has_alkylating_agent = False
                for reactant in reactants_smiles:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and (
                        mol.HasSubstructMatch(
                            Chem.MolFromSmarts("c1ccccc1C[Cl,Br,I,OH]")
                        )
                        or mol.HasSubstructMatch(Chem.MolFromSmarts("[C][Cl,Br,I]"))
                    ):
                        has_alkylating_agent = True
                        break

                # Check for tertiary amine in product
                product_mol = Chem.MolFromSmiles(product_smiles)
                has_tertiary_amine = False
                if product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[N;H0;!$(N=*)]([C])[C]")
                ):
                    has_tertiary_amine = True

                if has_secondary_amine and has_tertiary_amine:
                    print(f"Detected N-alkylation at depth {depth}")
                    has_n_alkylation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return has_n_alkylation

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
    This function detects a synthetic strategy involving late-stage fragment coupling
    via nucleophilic aromatic substitution.
    """
    late_stage_nas = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_nas

        if node["type"] == "reaction" and depth <= 1:  # Late stage (low depth)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nucleophilic aromatic substitution
                # Look for chloro-aromatic and amine reactants
                chloro_aromatic = False
                amine = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        # Check for chloro-aromatic
                        chloro_aromatic_pattern = Chem.MolFromSmarts("[c][Cl]")
                        if reactant_mol.HasSubstructMatch(chloro_aromatic_pattern):
                            chloro_aromatic = True

                        # Check for amine
                        amine_pattern = Chem.MolFromSmarts("[NH2]")
                        if reactant_mol.HasSubstructMatch(amine_pattern):
                            amine = True

                if chloro_aromatic and amine:
                    # Check if product has new C-N bond
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        c_n_bond_pattern = Chem.MolFromSmarts("[c][N]")
                        if product_mol.HasSubstructMatch(c_n_bond_pattern):
                            late_stage_nas = True
                            print(
                                "Detected late-stage nucleophilic aromatic substitution"
                            )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_nas

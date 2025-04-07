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
    Detects a reductive amination in the late stage of synthesis.
    """
    late_reductive_amination = False

    def dfs_traverse(node, depth=0):
        nonlocal late_reductive_amination

        if node["type"] == "reaction" and depth < 3:  # Late stage (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if (
                        all(mol is not None for mol in reactant_mols)
                        and product_mol is not None
                    ):
                        # Check for aldehyde and amine in reactants
                        aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
                        amine_pattern = Chem.MolFromSmarts("[NH]")

                        # Check for new C-N bond in product
                        cn_bond_pattern = Chem.MolFromSmarts("[CH2]-[N]")

                        has_aldehyde = any(
                            mol.HasSubstructMatch(aldehyde_pattern)
                            for mol in reactant_mols
                        )
                        has_amine = any(
                            mol.HasSubstructMatch(amine_pattern)
                            for mol in reactant_mols
                        )
                        has_cn_bond = product_mol.HasSubstructMatch(cn_bond_pattern)

                        if has_aldehyde and has_amine and has_cn_bond:
                            late_reductive_amination = True
                            print(
                                f"Late-stage reductive amination detected at depth {depth}"
                            )
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return late_reductive_amination

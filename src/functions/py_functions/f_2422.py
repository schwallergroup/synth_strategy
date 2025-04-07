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
    Detects a nitro reduction followed by cyclization to form a ring.
    """
    nitro_reduction_cyclization_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_cyclization_found

        if node["type"] == "reaction":
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
                        # Check for nitro group in reactants
                        nitro_pattern = Chem.MolFromSmarts("[c]-[N+](=O)[O-]")
                        lactam_pattern = Chem.MolFromSmarts("C1NC(=O)Cc2ccccc12")

                        has_nitro = any(
                            mol.HasSubstructMatch(nitro_pattern)
                            for mol in reactant_mols
                        )
                        has_lactam = product_mol.HasSubstructMatch(lactam_pattern)

                        # Check for ring formation
                        reactant_ring_count = sum(
                            Chem.GetSSSR(mol) for mol in reactant_mols
                        )
                        product_ring_count = Chem.GetSSSR(product_mol)

                        if (
                            has_nitro
                            and has_lactam
                            and product_ring_count > reactant_ring_count
                        ):
                            nitro_reduction_cyclization_found = True
                            print(
                                f"Nitro reduction with cyclization detected at depth {depth}"
                            )
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return nitro_reduction_cyclization_found

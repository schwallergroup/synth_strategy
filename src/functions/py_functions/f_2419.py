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
    Detects an O-alkylation in early stage followed by ring formation.
    """
    o_alkylation_depth = -1
    ring_formation_after_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal o_alkylation_depth, ring_formation_after_alkylation

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
                        # Check for O-alkylation (Williamson ether synthesis)
                        oh_pattern = Chem.MolFromSmarts("[OH]-[c]")
                        br_pattern = Chem.MolFromSmarts("Br[CH2][C](=O)[O]")
                        ether_pattern = Chem.MolFromSmarts("[c]-[O]-[CH2][C](=O)[O]")

                        has_oh = any(
                            mol.HasSubstructMatch(oh_pattern) for mol in reactant_mols
                        )
                        has_br = any(
                            mol.HasSubstructMatch(br_pattern) for mol in reactant_mols
                        )
                        has_ether = product_mol.HasSubstructMatch(ether_pattern)

                        if (
                            has_oh and has_br and has_ether and depth > 3
                        ):  # Early stage (high depth)
                            o_alkylation_depth = depth
                            print(f"O-alkylation detected at depth {depth}")

                        # Check for ring formation
                        reactant_ring_count = sum(
                            Chem.GetSSSR(mol) for mol in reactant_mols
                        )
                        product_ring_count = Chem.GetSSSR(product_mol)

                        if (
                            product_ring_count > reactant_ring_count
                            and o_alkylation_depth > -1
                            and depth < o_alkylation_depth
                        ):
                            ring_formation_after_alkylation = True
                            print(
                                f"Ring formation after O-alkylation detected at depth {depth}"
                            )
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return o_alkylation_depth > -1 and ring_formation_after_alkylation

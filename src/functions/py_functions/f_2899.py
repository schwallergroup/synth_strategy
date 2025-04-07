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
    Detects if the synthetic route involves an amide coupling reaction.
    """
    amide_coupling_found = False

    def dfs_traverse(node):
        nonlocal amide_coupling_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amide coupling
                if len(reactants) >= 2:  # Need at least two reactants
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    product_mol = Chem.MolFromSmiles(product) if product else None

                    if (
                        all(m is not None for m in reactant_mols)
                        and product_mol is not None
                    ):
                        # Check for carboxylic acid
                        carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
                        # Check for amine
                        amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
                        # Check for amide in product
                        amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3;H2,H1,H0]")

                        has_acid = any(
                            len(m.GetSubstructMatches(carboxylic_acid_pattern)) > 0
                            for m in reactant_mols
                        )
                        has_amine = any(
                            len(m.GetSubstructMatches(amine_pattern)) > 0
                            for m in reactant_mols
                        )
                        has_amide = product_mol.GetSubstructMatches(amide_pattern)

                        if has_acid and has_amine and has_amide:
                            amide_coupling_found = True
                            print("Found amide coupling reaction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return amide_coupling_found

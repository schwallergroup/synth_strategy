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
    Detects if the synthetic route uses trifluoroacetyl protection and deprotection of an amine.
    """
    protection_found = False
    deprotection_found = False

    def dfs_traverse(node):
        nonlocal protection_found, deprotection_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for protection: NH2 -> NH-C(=O)CF3
                if not protection_found:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    product_mol = Chem.MolFromSmiles(product) if product else None

                    if (
                        all(m is not None for m in reactant_mols)
                        and product_mol is not None
                    ):
                        # Check if any reactant has NH2 group
                        nh2_pattern = Chem.MolFromSmarts("[NX3H2]")
                        # Check if product has trifluoroacetyl group
                        tfacetyl_pattern = Chem.MolFromSmarts(
                            "[NX3]-C(=O)C([F])([F])([F])"
                        )

                        has_nh2 = any(
                            len(m.GetSubstructMatches(nh2_pattern)) > 0
                            for m in reactant_mols
                        )
                        has_tfacetyl = product_mol.GetSubstructMatches(tfacetyl_pattern)

                        if has_nh2 and has_tfacetyl:
                            protection_found = True
                            print("Found trifluoroacetyl protection reaction")

                # Check for deprotection: NH-C(=O)CF3 -> NH2
                if not deprotection_found:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    product_mol = Chem.MolFromSmiles(product) if product else None

                    if (
                        all(m is not None for m in reactant_mols)
                        and product_mol is not None
                    ):
                        # Check if any reactant has trifluoroacetyl group
                        tfacetyl_pattern = Chem.MolFromSmarts(
                            "[NX3]-C(=O)C([F])([F])([F])"
                        )
                        # Check if product has NH2 group
                        nh2_pattern = Chem.MolFromSmarts("[NX3H2]")

                        has_tfacetyl = any(
                            len(m.GetSubstructMatches(tfacetyl_pattern)) > 0
                            for m in reactant_mols
                        )
                        has_nh2 = product_mol.GetSubstructMatches(nh2_pattern)

                        if has_tfacetyl and has_nh2:
                            deprotection_found = True
                            print("Found trifluoroacetyl deprotection reaction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both protection and deprotection were found
    return protection_found and deprotection_found

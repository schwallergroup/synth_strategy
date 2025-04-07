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
    This function detects a synthetic strategy involving late-stage acylation
    of an amine on a heterocyclic scaffold.
    """
    acylation_found = False

    def dfs_traverse(node):
        nonlocal acylation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is the first (depth 0) reaction
            if any(
                child["type"] == "mol" and not child.get("in_stock", False)
                for child in node.get("children", [])
            ):
                # Check for acylation pattern: acyl chloride + amine → amide
                acyl_chloride_pattern = Chem.MolFromSmarts("[C](=[O])[Cl]")
                amine_pattern = Chem.MolFromSmarts("[NH2][c]")
                amide_pattern = Chem.MolFromSmarts("[NH]([c])[C](=[O])[#6]")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                if (
                    product_mol
                    and any(
                        mol and mol.HasSubstructMatch(acyl_chloride_pattern)
                        for mol in reactant_mols
                    )
                    and any(mol and mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols)
                    and product_mol.HasSubstructMatch(amide_pattern)
                ):
                    print("Late-stage acylation detected in first reaction")
                    acylation_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return acylation_found

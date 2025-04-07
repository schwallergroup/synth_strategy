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
    This function detects a synthetic strategy where a halogenated aromatic group
    is introduced late in the synthesis via reductive amination.
    """
    # Track if we found the pattern
    found_pattern = False

    def dfs_traverse(node):
        nonlocal found_pattern

        if node["type"] == "reaction" and not found_pattern:
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                if product_mol and all(reactant_mols):
                    # Check for reductive amination with halogenated aromatic
                    aldehyde_pattern = Chem.MolFromSmarts("[C;H1]=O")
                    sec_amine_pattern = Chem.MolFromSmarts("[N;H1]")
                    tert_amine_pattern = Chem.MolFromSmarts("[N;H0]")
                    halogen_pattern = Chem.MolFromSmarts("[c]~[Cl,Br,F,I]")

                    # Check if this is a reductive amination
                    is_reductive_amination = (
                        any(
                            r.HasSubstructMatch(aldehyde_pattern) for r in reactant_mols
                        )
                        and any(
                            r.HasSubstructMatch(sec_amine_pattern)
                            for r in reactant_mols
                        )
                        and product_mol.HasSubstructMatch(tert_amine_pattern)
                    )

                    # Check if halogenated aromatic is introduced
                    halogen_introduced = product_mol.HasSubstructMatch(
                        halogen_pattern
                    ) and any(
                        r.HasSubstructMatch(halogen_pattern) for r in reactant_mols
                    )

                    # Check if this is a late-stage reaction (depth 0 or 1)
                    if is_reductive_amination and halogen_introduced:
                        found_pattern = True
                        print(
                            "Detected late-stage reductive amination with halogenated aromatic"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    return found_pattern

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
    This function detects a sequence of functional group interconversions:
    ester → acid → acid chloride.
    """
    # Track the sequence of transformations
    transformations = []

    # SMARTS patterns
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    acid_chloride_pattern = Chem.MolFromSmarts("C(=O)Cl")

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                    product_mol = Chem.MolFromSmiles(product)

                    # Check for ester to acid conversion
                    has_ester_reactant = any(
                        mol is not None and mol.HasSubstructMatch(ester_pattern)
                        for mol in reactant_mols
                    )
                    has_acid_product = product_mol is not None and product_mol.HasSubstructMatch(
                        acid_pattern
                    )

                    if has_ester_reactant and has_acid_product:
                        transformations.append(("ester_to_acid", depth))
                        print(f"Found ester to acid conversion at depth {depth}")

                    # Check for acid to acid chloride conversion
                    has_acid_reactant = any(
                        mol is not None and mol.HasSubstructMatch(acid_pattern)
                        for mol in reactant_mols
                    )
                    has_acid_chloride_product = (
                        product_mol is not None
                        and product_mol.HasSubstructMatch(acid_chloride_pattern)
                    )

                    if has_acid_reactant and has_acid_chloride_product:
                        transformations.append(("acid_to_acid_chloride", depth))
                        print(f"Found acid to acid chloride conversion at depth {depth}")
                except:
                    print("Error processing SMILES in reaction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have the complete sequence
    has_ester_to_acid = any(t[0] == "ester_to_acid" for t in transformations)
    has_acid_to_acid_chloride = any(t[0] == "acid_to_acid_chloride" for t in transformations)

    strategy_present = has_ester_to_acid and has_acid_to_acid_chloride
    print(f"Functional group interconversion sequence: {strategy_present}")
    return strategy_present

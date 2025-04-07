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
    Detects convergent synthesis with late-stage C-N bond formation between
    a halogenated heterocycle and an amine-containing fragment.
    """
    found_late_stage_cn_coupling = False

    def dfs_traverse(node):
        nonlocal found_late_stage_cn_coupling

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            depth = node.get("metadata", {}).get("depth", -1)

            # Check if this is a late-stage reaction (depth 0 or 1)
            if depth <= 1:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                if len(reactants_smiles) >= 2:  # At least two reactants (convergent)
                    try:
                        reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                        product = Chem.MolFromSmiles(product_smiles)

                        # Check if any reactant has a halogen
                        has_halogen = any(
                            mol.GetSubstructMatches(Chem.MolFromSmarts("[F,Cl,Br,I]"))
                            for mol in reactants
                            if mol
                        )

                        # Check if any reactant has an amine
                        has_amine = any(
                            mol.GetSubstructMatches(Chem.MolFromSmarts("[NH2,NH]"))
                            for mol in reactants
                            if mol
                        )

                        # Check if product has a new C-N bond
                        if has_halogen and has_amine and product:
                            print(f"Found potential late-stage C-N coupling at depth {depth}")
                            found_late_stage_cn_coupling = True
                    except:
                        print("Error processing reaction SMILES")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_late_stage_cn_coupling

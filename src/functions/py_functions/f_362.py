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
    This function detects if the synthesis involves the formation of a piperazine ring
    via nucleophilic aromatic substitution.
    """
    piperazine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6][#6]1")
    chloro_aromatic_pattern = Chem.MolFromSmarts("[Cl]-[c]")

    has_piperazine_formation = False

    def dfs_traverse(node):
        nonlocal has_piperazine_formation

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                # Check if product has piperazine
                has_piperazine_in_product = (
                    product is not None
                    and product.HasSubstructMatch(piperazine_pattern)
                )

                # Check if any reactant has chloro-aromatic
                has_chloro_aromatic_in_reactants = any(
                    r is not None and r.HasSubstructMatch(chloro_aromatic_pattern)
                    for r in reactants
                )

                # Check if any reactant has piperazine precursor (secondary amine)
                has_secondary_amine_in_reactants = any(
                    r is not None and r.HasSubstructMatch(Chem.MolFromSmarts("[NH]"))
                    for r in reactants
                )

                # Check if this is a piperazine formation via nucleophilic aromatic substitution
                if (
                    has_piperazine_in_product
                    and has_chloro_aromatic_in_reactants
                    and has_secondary_amine_in_reactants
                ):
                    print(
                        "Found piperazine formation via nucleophilic aromatic substitution"
                    )
                    has_piperazine_formation = True
            except:
                print(f"Error processing reaction SMILES: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_piperazine_formation

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
    Detects if the synthesis route involves a nitration followed by reduction to amine
    followed by halogenation sequence.
    """
    # Track the sequence of transformations
    transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = (
                    Chem.MolFromSmiles(product_smiles) if product_smiles else None
                )

                if product_mol and reactant_mols:
                    # Check for nitration: product has nitro group that reactants don't have
                    nitro_pattern = Chem.MolFromSmarts("[#6]-[#7+](=[#8])[#8-]")
                    if product_mol.HasSubstructMatch(nitro_pattern) and not any(
                        r.HasSubstructMatch(nitro_pattern) for r in reactant_mols if r
                    ):
                        transformations.append(("nitration", depth))
                        print(f"Found nitration at depth {depth}")

                    # Check for reduction: reactant has nitro group, product has amine
                    amine_pattern = Chem.MolFromSmarts("[c]-[#7H2]")
                    if any(
                        r.HasSubstructMatch(nitro_pattern) for r in reactant_mols if r
                    ) and product_mol.HasSubstructMatch(amine_pattern):
                        transformations.append(("reduction", depth))
                        print(f"Found nitro reduction at depth {depth}")

                    # Check for halogenation of amine: reactant has amine, product has halide
                    bromo_pattern = Chem.MolFromSmarts("[c]-[#35]")
                    if any(
                        r.HasSubstructMatch(amine_pattern) for r in reactant_mols if r
                    ) and product_mol.HasSubstructMatch(bromo_pattern):
                        transformations.append(("halogenation", depth))
                        print(f"Found halogenation at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if the sequence exists in the correct order
    # Note: In retrosynthetic direction, the sequence would be halogenation -> reduction -> nitration
    has_nitration = any(t[0] == "nitration" for t in transformations)
    has_reduction = any(t[0] == "reduction" for t in transformations)
    has_halogenation = any(t[0] == "halogenation" for t in transformations)

    # Get the depths to check order
    nitration_depth = next((t[1] for t in transformations if t[0] == "nitration"), -1)
    reduction_depth = next((t[1] for t in transformations if t[0] == "reduction"), -1)
    halogenation_depth = next(
        (t[1] for t in transformations if t[0] == "halogenation"), -1
    )

    # Check if all transformations exist and are in the correct order (higher depth = earlier in synthesis)
    sequence_present = (
        has_nitration
        and has_reduction
        and has_halogenation
        and nitration_depth > reduction_depth > halogenation_depth
    )

    if sequence_present:
        print("Found complete nitration-reduction-halogenation sequence")

    return sequence_present

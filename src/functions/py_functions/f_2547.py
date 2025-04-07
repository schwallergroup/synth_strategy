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
    Detects if the synthesis follows a specific functional group transformation sequence:
    nitration → cyanation → reduction → hydrolysis → amidation → cyclization
    """
    # Track transformations in order of occurrence (higher depth = earlier in synthesis)
    transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactant_mols = [
                        Chem.MolFromSmiles(r)
                        for r in reactants_smiles
                        if Chem.MolFromSmiles(r)
                    ]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if product_mol and reactant_mols:
                        # Patterns for different transformations
                        nitro_pattern = Chem.MolFromSmarts("c[N+](=[O])[O-]")
                        cyano_pattern = Chem.MolFromSmarts("cC#N")
                        amine_pattern = Chem.MolFromSmarts("c[NH2]")
                        amide_pattern = Chem.MolFromSmarts("C(=[O])[NH2]")
                        acyl_pattern = Chem.MolFromSmarts("C(=[O])[NH]C")

                        # Check for nitration
                        if product_mol.HasSubstructMatch(nitro_pattern) and not any(
                            r.HasSubstructMatch(nitro_pattern) for r in reactant_mols
                        ):
                            transformations.append(("nitration", depth))
                            print(f"Nitration detected at depth {depth}")

                        # Check for cyanation
                        if product_mol.HasSubstructMatch(cyano_pattern) and not any(
                            r.HasSubstructMatch(cyano_pattern) for r in reactant_mols
                        ):
                            transformations.append(("cyanation", depth))
                            print(f"Cyanation detected at depth {depth}")

                        # Check for reduction (nitro to amine)
                        if product_mol.HasSubstructMatch(amine_pattern) and any(
                            r.HasSubstructMatch(nitro_pattern) for r in reactant_mols
                        ):
                            transformations.append(("reduction", depth))
                            print(f"Reduction detected at depth {depth}")

                        # Check for hydrolysis (nitrile to amide)
                        if product_mol.HasSubstructMatch(amide_pattern) and any(
                            r.HasSubstructMatch(cyano_pattern) for r in reactant_mols
                        ):
                            transformations.append(("hydrolysis", depth))
                            print(f"Hydrolysis detected at depth {depth}")

                        # Check for amidation
                        if product_mol.HasSubstructMatch(acyl_pattern) and not any(
                            r.HasSubstructMatch(acyl_pattern) for r in reactant_mols
                        ):
                            transformations.append(("amidation", depth))
                            print(f"Amidation detected at depth {depth}")

                        # Check for cyclization (ring count increase)
                        reactant_ring_count = sum(
                            Chem.GetSSSR(r) for r in reactant_mols
                        )
                        product_ring_count = Chem.GetSSSR(product_mol)
                        if product_ring_count > reactant_ring_count:
                            transformations.append(("cyclization", depth))
                            print(f"Cyclization detected at depth {depth}")

                except Exception as e:
                    print(f"Error in transformation detection: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Sort transformations by depth (higher depth = earlier in synthesis)
    transformations.sort(key=lambda x: x[1], reverse=True)

    # Extract just the transformation types in sequence
    sequence = [t[0] for t in transformations]
    print(f"Detected transformation sequence: {sequence}")

    # Check if the sequence matches our target pattern
    # We're looking for these transformations in order (some might be missing)
    target_sequence = [
        "nitration",
        "cyanation",
        "reduction",
        "hydrolysis",
        "amidation",
        "cyclization",
    ]

    # Check if sequence follows the target order (allowing for missing steps)
    if len(sequence) >= 3:  # Require at least 3 of the steps
        last_index = -1
        matches = 0
        for step in sequence:
            if step in target_sequence:
                current_index = target_sequence.index(step)
                if current_index > last_index:
                    matches += 1
                    last_index = current_index

        # Return True if at least 3 steps match and are in the correct order
        return matches >= 3

    return False

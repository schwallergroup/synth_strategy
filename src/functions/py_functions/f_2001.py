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
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    Detects a linear fragment assembly strategy where fragments are added sequentially
    without convergent steps, building up to a final cyclization.
    """
    # Track the pattern
    is_linear = True
    has_final_cyclization = False
    step_count = 0
    first_step_processed = False

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, has_final_cyclization, step_count, first_step_processed

        if node["type"] == "reaction":
            step_count += 1

            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Processing reaction at depth {depth}: {rsmi}")

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if not reactant_mols or not product_mol:
                print("Could not parse molecules in reaction")
                return

            # Check if this is a convergent step (more than 2 significant fragments combining)
            significant_fragments = 0
            for mol in reactant_mols:
                # Consider only non-trivial fragments (more than 5 atoms)
                if mol.GetNumAtoms() > 5:
                    significant_fragments += 1

            print(
                f"Found {significant_fragments} significant fragments at depth {depth}"
            )

            # In linear assembly, we should have at most 2 significant fragments
            # (the growing chain and the new fragment being added)
            if significant_fragments > 2:
                is_linear = False
                print(
                    f"Found convergent step at depth {depth} with {significant_fragments} significant fragments"
                )

            # Check for cyclization in final step (first reaction encountered in retrosynthetic traversal)
            if not first_step_processed:
                first_step_processed = True

                # Count rings in reactants and product
                reactant_ring_count = sum(
                    len(Chem.GetSSSR(mol)) for mol in reactant_mols
                )
                product_ring_count = len(Chem.GetSSSR(product_mol))

                print(
                    f"Ring count: reactants={reactant_ring_count}, product={product_ring_count}"
                )

                # Check if this is a cyclization reaction
                if product_ring_count > reactant_ring_count:
                    has_final_cyclization = True
                    print("Found cyclization in final step (ring count increased)")
                # Also check for known cyclization reaction types
                elif checker.check_reaction("Formation of NOS Heterocycles", rsmi):
                    has_final_cyclization = True
                    print("Found cyclization in final step (NOS Heterocycle formation)")
                elif any(
                    checker.check_reaction(rxn, rsmi)
                    for rxn in [
                        "Paal-Knorr pyrrole synthesis",
                        "Intramolecular transesterification/Lactone formation",
                        "Intramolecular amination (heterocycle formation)",
                    ]
                ):
                    has_final_cyclization = True
                    print(
                        "Found cyclization in final step (known cyclization reaction)"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Need at least 3 steps to be considered a meaningful linear strategy
    is_meaningful_linear = is_linear and step_count >= 3

    print(
        f"Analysis complete: linear={is_linear}, steps={step_count}, final_cyclization={has_final_cyclization}"
    )

    if is_meaningful_linear and has_final_cyclization:
        print("Detected linear fragment assembly with final cyclization strategy")
        return True
    elif is_meaningful_linear:
        print("Detected linear fragment assembly but no final cyclization")
        return False
    else:
        print("Not a linear fragment assembly strategy")
        return False

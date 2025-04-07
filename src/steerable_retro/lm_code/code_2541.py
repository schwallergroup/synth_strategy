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
    Detects if the synthesis involves a late-stage heterocyclic ring formation
    in the final or penultimate step.
    """
    ring_formation_detected = False

    # List of common heterocyclic rings to check
    heterocyclic_rings = [
        "furan",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "indole",
        "benzofuran",
        "benzothiophene",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
        "thiophene",
        "morpholine",
        "piperidine",
        "piperazine",
        "pyrrolidine",
        "tetrahydrofuran",
        "dioxane",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal ring_formation_detected

        if node["type"] == "reaction" and depth <= 1:  # Final or penultimate step
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Convert SMILES to RDKit molecules
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if product_mol and all(r for r in reactant_mols):
                        # Count rings in reactants and product
                        reactant_ring_count = sum(len(Chem.GetSSSR(r)) for r in reactant_mols)
                        product_ring_count = len(Chem.GetSSSR(product_mol))

                        # Check if product has more rings than reactants combined
                        if product_ring_count > reactant_ring_count:
                            print(
                                f"Ring count increased: {reactant_ring_count} → {product_ring_count}"
                            )

                            # Check for heterocyclic rings in product that aren't in reactants
                            for ring_name in heterocyclic_rings:
                                if checker.check_ring(ring_name, product_smiles):
                                    # Verify this ring wasn't present in any reactant
                                    ring_is_new = True
                                    for reactant in reactants_smiles:
                                        if checker.check_ring(ring_name, reactant):
                                            ring_is_new = False
                                            break

                                    if ring_is_new:
                                        print(
                                            f"Late-stage heterocyclic ring formation detected: {ring_name} at depth {depth}"
                                        )
                                        ring_formation_detected = True
                                        break

                            # If no specific heterocycle was identified but rings increased,
                            # check for general heteroatoms in rings
                            if not ring_formation_detected:
                                # Check for any heterocyclic structure
                                hetero_pattern = Chem.MolFromSmarts(
                                    "[r;!#6]"
                                )  # Ring atom that's not carbon
                                if product_mol.HasSubstructMatch(hetero_pattern):
                                    # Check if this pattern is new
                                    hetero_is_new = True
                                    for r_mol in reactant_mols:
                                        if r_mol and r_mol.HasSubstructMatch(hetero_pattern):
                                            # Need to check if the specific heterocyclic system is new
                                            # This is a simplified check
                                            if len(
                                                r_mol.GetSubstructMatches(hetero_pattern)
                                            ) >= len(
                                                product_mol.GetSubstructMatches(hetero_pattern)
                                            ):
                                                hetero_is_new = False
                                                break

                                    if hetero_is_new:
                                        print(
                                            f"Generic heterocyclic ring formation detected at depth {depth}"
                                        )
                                        ring_formation_detected = True
                except Exception as e:
                    print(f"Error in ring formation detection: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return ring_formation_detected

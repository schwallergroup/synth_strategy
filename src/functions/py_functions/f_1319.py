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
    This function detects a late-stage cyclization strategy where a ring is formed
    in one of the final steps of the synthesis.
    """
    cyclization_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal cyclization_detected

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Focus on late-stage reactions (depth 0 or 1)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Process each reactant separately
                    reactants_list = reactants_smiles.split(".")
                    reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_list]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if (
                        all(mol is not None for mol in reactants_mols)
                        and product_mol is not None
                    ):
                        # Count rings in reactants and product
                        reactant_ring_count = sum(
                            len(mol.GetRingInfo().AtomRings()) for mol in reactants_mols
                        )
                        product_ring_count = len(product_mol.GetRingInfo().AtomRings())

                        print(
                            f"Depth: {depth}, Reactant rings: {reactant_ring_count}, Product rings: {product_ring_count}"
                        )

                        # Check if a new ring was formed
                        if product_ring_count > reactant_ring_count:
                            # Verify it's a recognized ring structure
                            ring_types = [
                                "furan",
                                "pyran",
                                "pyrrole",
                                "pyridine",
                                "benzene",
                                "naphthalene",
                                "cyclopentane",
                                "cyclohexane",
                                "indole",
                                "quinoline",
                                "thiophene",
                                "oxazole",
                                "thiazole",
                                "imidazole",
                                "pyrazole",
                                "triazole",
                                "tetrazole",
                                "piperidine",
                                "piperazine",
                                "morpholine",
                            ]

                            for ring_type in ring_types:
                                if checker.check_ring(
                                    ring_type, product_smiles
                                ) and not any(
                                    checker.check_ring(ring_type, r)
                                    for r in reactants_list
                                ):
                                    print(
                                        f"Late-stage cyclization detected at depth {depth}: {ring_type} ring formed"
                                    )
                                    cyclization_detected = True
                                    break

                            # If no specific ring identified but ring count increased, still mark as cyclization
                            if (
                                not cyclization_detected
                                and product_ring_count > reactant_ring_count
                            ):
                                print(
                                    f"Late-stage cyclization detected at depth {depth}: unspecified ring formed"
                                )
                                cyclization_detected = True
                except Exception as e:
                    print(f"Error processing SMILES: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return cyclization_detected

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
    Detects ring formation in the early stages of synthesis (high depth).
    """
    found_ring_formation = False

    # List of common ring types to check
    ring_types = [
        "furan",
        "pyran",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "cyclopropane",
        "cyclobutane",
        "cyclopentane",
        "cyclohexane",
        "benzene",
        "indole",
        "quinoline",
    ]

    # List of common ring-forming reaction types
    ring_forming_reactions = [
        "Diels-Alder",
        "Paal-Knorr pyrrole synthesis",
        "Pictet-Spengler",
        "Fischer indole",
        "Friedlaender chinoline",
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_ring_formation

        if node["type"] == "reaction" and depth >= 4:  # Early in synthesis
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if not product or not all(reactants):
                    return

                # Check if this is a known ring-forming reaction
                for rxn_type in ring_forming_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected early ring-forming reaction: {rxn_type} at depth {depth}")
                        found_ring_formation = True
                        return

                # Count rings in reactants and product
                reactant_rings = sum([r.GetRingInfo().NumRings() for r in reactants])
                product_rings = product.GetRingInfo().NumRings()

                if product_rings > reactant_rings:
                    print(f"Detected early ring formation at depth {depth}")
                    print(
                        f"Rings in reactants: {reactant_rings}, Rings in product: {product_rings}"
                    )

                    # Check which specific rings are formed
                    for ring_type in ring_types:
                        # Check if ring exists in product but not in any reactant
                        ring_in_product = checker.check_ring(ring_type, product_smiles)
                        ring_in_reactants = any(
                            checker.check_ring(ring_type, r) for r in reactants_smiles
                        )

                        if ring_in_product and not ring_in_reactants:
                            print(f"Formed {ring_type} ring at depth {depth}")
                            found_ring_formation = True
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Early ring formation detected: {found_ring_formation}")
    return found_ring_formation

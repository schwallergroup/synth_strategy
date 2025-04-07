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
    Detects a synthetic strategy involving a late-stage Suzuki coupling with preparatory borylation
    and earlier ring formation.
    """
    has_suzuki_coupling = False
    has_borylation = False
    has_ring_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_suzuki_coupling, has_borylation, has_ring_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Suzuki coupling (late stage, depth 0-1)
                if depth <= 1:
                    # Check directly for Suzuki coupling reaction
                    if (
                        checker.check_reaction(
                            "Suzuki coupling with boronic acids", rsmi
                        )
                        or checker.check_reaction(
                            "Suzuki coupling with boronic acids OTf", rsmi
                        )
                        or checker.check_reaction(
                            "Suzuki coupling with boronic esters", rsmi
                        )
                        or checker.check_reaction(
                            "Suzuki coupling with boronic esters OTf", rsmi
                        )
                    ):
                        has_suzuki_coupling = True
                        print(f"Found late-stage Suzuki coupling at depth {depth}")
                    # Fallback check using functional groups
                    elif any(
                        checker.check_fg("Boronic acid", r)
                        or checker.check_fg("Boronic ester", r)
                        for r in reactants
                    ) and any(
                        checker.check_fg("Aromatic halide", r) for r in reactants
                    ):
                        has_suzuki_coupling = True
                        print(
                            f"Found late-stage Suzuki coupling (FG check) at depth {depth}"
                        )

                # Check for borylation (mid stage, depth 1-2)
                elif depth <= 3:
                    # Check for preparation of boronic acid/ester
                    if (
                        checker.check_reaction("Preparation of boronic acids", rsmi)
                        or checker.check_reaction("Preparation of boronic ethers", rsmi)
                        or checker.check_reaction(
                            "Preparation of boronic acids from trifluoroborates", rsmi
                        )
                    ):
                        has_borylation = True
                        print(f"Found borylation step at depth {depth}")
                    # Fallback check using functional groups
                    elif (
                        checker.check_fg("Boronic acid", product)
                        or checker.check_fg("Boronic ester", product)
                    ) and any(
                        checker.check_fg("Aromatic halide", r) for r in reactants
                    ):
                        has_borylation = True
                        print(f"Found borylation step (FG check) at depth {depth}")

                # Check for ring formation (early stage, depth 2+)
                if depth >= 2:
                    # Check for common ring-forming reactions
                    ring_forming_reactions = [
                        "Paal-Knorr pyrrole synthesis",
                        "Formation of NOS Heterocycles",
                        "Benzothiazole formation",
                        "Benzoxazole formation",
                        "Benzimidazole formation",
                        "Diels-Alder",
                        "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                    ]

                    if any(
                        checker.check_reaction(rxn, rsmi)
                        for rxn in ring_forming_reactions
                    ):
                        has_ring_formation = True
                        print(f"Found ring formation reaction at depth {depth}")
                        return

                    # Check for common rings in product but not in all reactants
                    common_rings = [
                        "pyrrole",
                        "pyridine",
                        "furan",
                        "thiophene",
                        "benzene",
                        "indole",
                        "pyrazole",
                        "imidazole",
                        "oxazole",
                        "thiazole",
                        "pyrimidine",
                        "quinoline",
                        "isoquinoline",
                        "benzothiazole",
                        "benzoxazole",
                        "benzimidazole",
                    ]

                    # Check if product has rings that aren't in all reactants
                    product_rings = [
                        ring
                        for ring in common_rings
                        if checker.check_ring(ring, product)
                    ]
                    if product_rings:
                        # Check if any ring in product is not present in all reactants
                        for ring in product_rings:
                            if not all(
                                checker.check_ring(ring, r)
                                for r in reactants
                                if Chem.MolFromSmiles(r)
                            ):
                                has_ring_formation = True
                                print(f"Found ring formation at depth {depth}: {ring}")
                                break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if all key elements of the strategy are present
    result = has_suzuki_coupling and has_borylation and has_ring_formation
    print(
        f"Strategy detection result: Suzuki={has_suzuki_coupling}, Borylation={has_borylation}, Ring formation={has_ring_formation}"
    )
    return result

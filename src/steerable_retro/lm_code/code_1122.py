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
    Detects a strategy involving late-stage arylation via cross-coupling,
    particularly Suzuki-type coupling with boronic acids.
    """
    late_arylation_found = False
    depth_threshold = 2  # Consider "late stage" if depth <= 2

    def dfs_traverse(node, depth=0):
        nonlocal late_arylation_found

        if node["type"] == "reaction" and depth <= depth_threshold:
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is an arylation cross-coupling reaction
                is_arylation = (
                    checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                    or checker.check_reaction("Buchwald-Hartwig", rsmi)
                    or checker.check_reaction(
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                    )
                    or checker.check_reaction("Heck terminal vinyl", rsmi)
                    or checker.check_reaction("Negishi coupling", rsmi)
                    or checker.check_reaction("Stille reaction_aryl", rsmi)
                    or checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi)
                    or checker.check_reaction("Ullmann condensation", rsmi)
                    or checker.check_reaction("Ullmann-Goldberg Substitution aryl alcohol", rsmi)
                    or checker.check_reaction("Hiyama-Denmark Coupling", rsmi)
                    or checker.check_reaction("Kumada cross-coupling", rsmi)
                )

                print(f"Is arylation cross-coupling: {is_arylation}")

                # Check for coupling partners in reactants
                coupling_partner_present = any(
                    checker.check_fg("Boronic acid", r)
                    or checker.check_fg("Boronic ester", r)
                    or checker.check_fg("Aromatic halide", r)
                    or checker.check_fg("Triflate", r)
                    for r in reactants
                )

                print(f"Has coupling partner: {coupling_partner_present}")

                # Check for aromatic rings in product
                aromatic_rings = [
                    "benzene",
                    "naphthalene",
                    "anthracene",
                    "pyridine",
                    "quinoline",
                    "isoquinoline",
                    "indole",
                    "benzothiophene",
                    "benzoxazole",
                    "benzimidazole",
                    "furan",
                    "thiophene",
                    "pyrazole",
                    "imidazole",
                    "oxazole",
                    "thiazole",
                    "pyrimidine",
                ]

                has_aromatic_ring = any(
                    checker.check_ring(ring, product) for ring in aromatic_rings
                )

                print(f"Has aromatic ring in product: {has_aromatic_ring}")

                # Check if a new C-C or C-N or C-O bond is formed between aromatic rings
                # This is a heuristic check for arylation reactions that might not be captured by reaction types
                new_aryl_bond_formed = False

                # If we have aromatic halides in reactants and aromatic rings in product
                # it's likely an arylation reaction
                aromatic_halide_present = any(
                    checker.check_fg("Aromatic halide", r) for r in reactants
                )

                if aromatic_halide_present and has_aromatic_ring:
                    print("Aromatic halide present in reactants and aromatic ring in product")
                    new_aryl_bond_formed = True

                if (
                    is_arylation or coupling_partner_present or new_aryl_bond_formed
                ) and has_aromatic_ring:
                    print(f"Found late-stage arylation via cross-coupling at depth {depth}")
                    print(
                        f"Is arylation: {is_arylation}, Has coupling partner: {coupling_partner_present}"
                    )
                    late_arylation_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return late_arylation_found

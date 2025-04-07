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
    Detects a strategy involving late-stage hydrazide formation (in the final or penultimate step).
    """
    has_late_stage_hydrazide = False

    def dfs_traverse(node):
        nonlocal has_late_stage_hydrazide

        if node["type"] == "reaction":
            # Check if this is a late-stage step (depth 0 or 1)
            # Treat missing depth as 0 (latest stage)
            if node.get("metadata", {}).get("depth", 0) <= 1:
                try:
                    if "rsmi" not in node.get("metadata", {}):
                        return

                    rsmi = node["metadata"]["rsmi"]
                    reactants_part = rsmi.split(">")[0]
                    reactants = reactants_part.split(".")
                    product = rsmi.split(">")[-1]

                    # Check if product contains hydrazide or hydrazone amide
                    product_has_hydrazide = checker.check_fg(
                        "Acylhydrazine", product
                    ) or checker.check_fg("Hydrazone amide", product)

                    # Check if any reactant contains hydrazide or hydrazone amide
                    reactants_have_hydrazide = any(
                        checker.check_fg("Acylhydrazine", reactant)
                        or checker.check_fg("Hydrazone amide", reactant)
                        for reactant in reactants
                    )

                    # Check for hydrazine in reactants
                    has_hydrazine = any(
                        checker.check_fg("Hydrazine", reactant)
                        for reactant in reactants
                    )

                    # Check for acyl halide, ester, carboxylic acid, or anhydride in reactants
                    has_acyl_source = any(
                        checker.check_fg("Acyl halide", reactant)
                        or checker.check_fg("Ester", reactant)
                        or checker.check_fg("Carboxylic acid", reactant)
                        or checker.check_fg("Anhydride", reactant)
                        for reactant in reactants
                    )

                    # Verify hydrazide formation: product has hydrazide, reactants don't, and appropriate reactants present
                    is_hydrazide_formation = False

                    if (
                        product_has_hydrazide
                        and not reactants_have_hydrazide
                        and has_hydrazine
                        and has_acyl_source
                    ):
                        is_hydrazide_formation = True

                    # Check for specific reactions that form hydrazides
                    if (
                        checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                            rsmi,
                        )
                        or checker.check_reaction("Aminolysis of esters", rsmi)
                        or checker.check_reaction(
                            "Carboxylic acid with primary amine to amide", rsmi
                        )
                        or checker.check_reaction("Acylation of primary amines", rsmi)
                        or checker.check_reaction(
                            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                            rsmi,
                        )
                        or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                    ) and has_hydrazine:
                        is_hydrazide_formation = True

                    # Also check for direct conversion of nitrile to hydrazide
                    if (
                        checker.check_reaction(
                            "Nitrile and hydrogen peroxide to amide", rsmi
                        )
                        and has_hydrazine
                    ):
                        is_hydrazide_formation = True

                    # Check for amide formation reactions that might form hydrazides
                    if (
                        checker.check_reaction(
                            "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                            rsmi,
                        )
                        or checker.check_reaction(
                            "Carboxylic acid to amide conversion", rsmi
                        )
                    ) and has_hydrazine:
                        is_hydrazide_formation = True

                    if is_hydrazide_formation:
                        has_late_stage_hydrazide = True
                        print(
                            f"Found hydrazide formation in late-stage step (depth {node['metadata'].get('depth', 0)}): {rsmi}"
                        )

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Late-stage hydrazide formation strategy detected: {has_late_stage_hydrazide}"
    )
    return has_late_stage_hydrazide

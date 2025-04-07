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
    Detects if the synthesis uses a late-stage thioether coupling strategy
    to connect two complex fragments.
    """
    late_stage_thioether_coupling_found = False

    def dfs_traverse(node):
        nonlocal late_stage_thioether_coupling_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Extract depth from ID - default to late stage (0) if not specified
            depth = 0  # Default to low depth (late stage)
            id_str = node.get("metadata", {}).get("ID", "")
            if "Depth:" in id_str:
                try:
                    depth = int(id_str.split("Depth:")[1].strip().split()[0])
                except (ValueError, IndexError):
                    pass

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Consider depths 0, 1, and 2 as late-stage
            if depth <= 2:
                # Check if this is a thioether coupling reaction
                is_thioether_reaction = (
                    checker.check_reaction("thioether_nucl_sub", rsmi)
                    or checker.check_reaction("S-alkylation of thiols", rsmi)
                    or checker.check_reaction(
                        "S-alkylation of thiols with alcohols", rsmi
                    )
                    or checker.check_reaction("S-alkylation of thiols (ethyl)", rsmi)
                )

                # Check if product contains a thioether (monosulfide)
                has_thioether_product = checker.check_fg("Monosulfide", product)

                # Check if thioether is newly formed (not present in reactants)
                thioether_in_reactants = any(
                    checker.check_fg("Monosulfide", r) for r in reactants
                )

                # Check for thiol reactant (required for thioether formation)
                has_thiol_reactant = any(
                    checker.check_fg("Aromatic thiol", r)
                    or checker.check_fg("Aliphatic thiol", r)
                    for r in reactants
                )

                # Check for thiocarbonyl reactant (can be converted to thioether)
                has_thiocarbonyl_reactant = any(
                    checker.check_fg("Thiocarbonyl", r) for r in reactants
                )

                # Check if coupling connects complex fragments
                complex_fragments = 0
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if (
                        mol and mol.GetNumHeavyAtoms() >= 8
                    ):  # Define "complex" as having 8+ heavy atoms
                        complex_fragments += 1
                connects_complex_fragments = complex_fragments >= 1

                print(f"  Thioether reaction: {is_thioether_reaction}")
                print(f"  Has thioether product: {has_thioether_product}")
                print(f"  Thioether in reactants: {thioether_in_reactants}")
                print(f"  Has thiol reactant: {has_thiol_reactant}")
                print(f"  Has thiocarbonyl reactant: {has_thiocarbonyl_reactant}")
                print(f"  Complex fragments: {complex_fragments}")

                # Check for the specific case: thiocarbonyl to thioether conversion
                if (
                    has_thioether_product
                    and has_thiocarbonyl_reactant
                    and not thioether_in_reactants
                    and connects_complex_fragments
                ):
                    print(
                        f"  Special case: Thiocarbonyl to thioether conversion at depth {depth}"
                    )
                    print(f"  Product: {product}")
                    print(f"  Reactants: {reactants}")
                    late_stage_thioether_coupling_found = True

                # If all conditions are met, we've found a late-stage thioether coupling
                elif (
                    is_thioether_reaction
                    and has_thioether_product
                    and (has_thiol_reactant or has_thiocarbonyl_reactant)
                    and not thioether_in_reactants
                    and connects_complex_fragments
                ):
                    print(f"  Late-stage thioether coupling found at depth {depth}")
                    late_stage_thioether_coupling_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(
        f"Late-stage thioether coupling detected: {late_stage_thioether_coupling_found}"
    )
    return late_stage_thioether_coupling_found

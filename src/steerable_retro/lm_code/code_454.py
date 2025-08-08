#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
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
    This function detects a synthetic strategy involving nitration-reduction sequence
    followed by a late-stage amide coupling.
    """
    # Track if we found the key elements of the strategy
    found_amide_coupling = False
    found_nitro_reduction = False
    amide_coupling_depth = float("inf")
    nitro_reduction_depth = float("inf")

    # Track molecules with reduced nitro groups
    reduced_nitro_molecules = set()

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_coupling, found_nitro_reduction, amide_coupling_depth, nitro_reduction_depth

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for amide coupling reactions at any depth
            if (
                checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                )
                or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                or checker.check_reaction(
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                )
                or checker.check_reaction("Ester with primary amine to amide", rsmi)
                or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                or checker.check_reaction("Schotten-Baumann_amide", rsmi)
            ):
                found_amide_coupling = True
                amide_coupling_depth = min(amide_coupling_depth, depth)
                print(f"Found amide coupling at depth {depth}")
            else:
                # Fallback check for amide coupling
                has_carboxylic_acid = any(
                    checker.check_fg("Carboxylic acid", r) for r in reactants_smiles if r
                )
                has_acyl_halide = any(
                    checker.check_fg("Acyl halide", r) for r in reactants_smiles if r
                )
                has_ester = any(checker.check_fg("Ester", r) for r in reactants_smiles if r)
                has_amine = any(
                    checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                    for r in reactants_smiles
                    if r
                )
                has_amide_product = (
                    checker.check_fg("Primary amide", product_smiles)
                    or checker.check_fg("Secondary amide", product_smiles)
                    or checker.check_fg("Tertiary amide", product_smiles)
                )

                if (
                    has_amide_product
                    and (has_carboxylic_acid or has_acyl_halide or has_ester)
                    and has_amine
                ):
                    found_amide_coupling = True
                    amide_coupling_depth = min(amide_coupling_depth, depth)
                    print(f"Found amide coupling (fallback check) at depth {depth}")

            # Check for nitro reduction
            if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                found_nitro_reduction = True
                nitro_reduction_depth = min(nitro_reduction_depth, depth)
                print(f"Found nitro reduction at depth {depth}")
                # Add the product to our set of molecules with reduced nitro groups
                reduced_nitro_molecules.add(product_smiles)
            else:
                # Fallback check if reaction type not directly available
                # Check if reactant has nitro group and product has amine
                has_nitro_reactant = any(
                    checker.check_fg("Nitro group", r) for r in reactants_smiles if r
                )
                has_amine_product = (
                    (
                        checker.check_fg("Primary amine", product_smiles)
                        or checker.check_fg("Secondary amine", product_smiles)
                    )
                    if product_smiles
                    else False
                )

                # Make sure the nitro group is actually being reduced to an amine
                # by checking that the nitro group is not present in the product
                if (
                    has_nitro_reactant
                    and has_amine_product
                    and not checker.check_fg("Nitro group", product_smiles)
                ):
                    found_nitro_reduction = True
                    nitro_reduction_depth = min(nitro_reduction_depth, depth)
                    print(f"Found nitro reduction (fallback check) at depth {depth}")
                    # Add the product to our set of molecules with reduced nitro groups
                    reduced_nitro_molecules.add(product_smiles)

                # Also check for nitration reactions (adding a nitro group)
                if (
                    checker.check_reaction("Aromatic nitration with HNO3", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO3 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO2 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with alkyl NO2", rsmi)
                    or checker.check_reaction("Non-aromatic nitration with HNO3", rsmi)
                ):
                    print(f"Found nitration reaction at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the strategy is present: both reactions found and amide coupling occurs after nitro reduction
    # For late-stage amide coupling, the amide coupling should be at a lower depth (closer to final product)
    # than the nitro reduction
    strategy_present = (
        found_amide_coupling
        and found_nitro_reduction
        and amide_coupling_depth < nitro_reduction_depth
    )

    if strategy_present:
        print(
            f"Detected late-stage amide coupling (depth {amide_coupling_depth}) with preceding nitro reduction (depth {nitro_reduction_depth})"
        )
    else:
        if found_amide_coupling:
            print(f"Found amide coupling at depth {amide_coupling_depth}")
        if found_nitro_reduction:
            print(f"Found nitro reduction at depth {nitro_reduction_depth}")
        if found_amide_coupling and found_nitro_reduction:
            print(
                f"Both reactions found but in wrong order: amide at depth {amide_coupling_depth}, nitro reduction at depth {nitro_reduction_depth}"
            )
        else:
            print("Did not find both required reactions")

    return strategy_present

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
    This function detects a strategy where a nitro group is reduced to an amine
    before a coupling reaction occurs later in the synthesis.
    """
    # Track molecules with nitro groups and their reduced amine products
    molecule_info = {}
    # Track reactions and their depths
    nitro_reductions = {}  # {product_smiles: depth}
    coupling_reactions = {}  # {product_smiles: (depth, reactant_smiles)}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check if molecule has nitro group
            if checker.check_fg("Nitro group", mol_smiles):
                molecule_info[mol_smiles] = molecule_info.get(mol_smiles, {})
                molecule_info[mol_smiles]["has_nitro"] = True
                molecule_info[mol_smiles]["depth"] = depth
                print(f"Found molecule with nitro group at depth {depth}: {mol_smiles}")

            # Check if molecule has primary amine
            if checker.check_fg("Primary amine", mol_smiles):
                molecule_info[mol_smiles] = molecule_info.get(mol_smiles, {})
                molecule_info[mol_smiles]["has_amine"] = True
                molecule_info[mol_smiles]["depth"] = depth
                print(
                    f"Found molecule with primary amine at depth {depth}: {mol_smiles}"
                )

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro reduction
                try:
                    reactant_has_nitro = any(
                        checker.check_fg("Nitro group", reactant)
                        for reactant in reactants
                    )
                    product_has_amine = checker.check_fg("Primary amine", product)

                    # Check if this is a nitro reduction reaction
                    is_nitro_reduction = checker.check_reaction(
                        "Reduction of nitro groups to amines", rsmi
                    )

                    # If not explicitly a nitro reduction, check if nitro disappears and amine appears
                    if (
                        not is_nitro_reduction
                        and reactant_has_nitro
                        and product_has_amine
                    ):
                        # Check if nitro group is gone in product
                        if not checker.check_fg("Nitro group", product):
                            is_nitro_reduction = True

                    if reactant_has_nitro and product_has_amine and is_nitro_reduction:
                        print(f"Found nitro reduction at depth {depth}, rsmi: {rsmi}")
                        nitro_reductions[product] = depth

                        # Update molecule info
                        for reactant in reactants:
                            if checker.check_fg("Nitro group", reactant):
                                molecule_info[reactant] = molecule_info.get(
                                    reactant, {}
                                )
                                molecule_info[reactant]["has_nitro"] = True
                                molecule_info[reactant]["reduced_to"] = product
                                molecule_info[reactant]["depth"] = depth

                        molecule_info[product] = molecule_info.get(product, {})
                        molecule_info[product]["has_amine"] = True
                        molecule_info[product]["reduced_from"] = [
                            r for r in reactants if checker.check_fg("Nitro group", r)
                        ]
                        molecule_info[product]["depth"] = depth
                except Exception as e:
                    print(f"Error checking nitro reduction: {e}")

                # Check for coupling reactions involving amines
                try:
                    # Check for various coupling reactions
                    coupling_reactions_list = [
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                        "Suzuki coupling with boronic acids",
                        "Sonogashira alkyne_aryl halide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Carboxylic acid with primary amine to amide",
                        "Ester with primary amine to amide",
                    ]

                    is_coupling = any(
                        checker.check_reaction(rxn_type, rsmi)
                        for rxn_type in coupling_reactions_list
                    )

                    # Also check for general coupling patterns if specific reaction types aren't detected
                    if not is_coupling:
                        # Check for amide formation (common coupling product)
                        if checker.check_fg(
                            "Primary amide", product
                        ) or checker.check_fg("Secondary amide", product):
                            reactant_has_amine = any(
                                checker.check_fg("Primary amine", reactant)
                                for reactant in reactants
                            )
                            reactant_has_acid = any(
                                checker.check_fg("Carboxylic acid", reactant)
                                for reactant in reactants
                            )
                            if reactant_has_amine and reactant_has_acid:
                                is_coupling = True

                    if is_coupling:
                        # Check if any reactant has primary amine
                        amine_reactants = [
                            r for r in reactants if checker.check_fg("Primary amine", r)
                        ]

                        if amine_reactants:
                            print(
                                f"Found coupling reaction using amine at depth {depth}, rsmi: {rsmi}"
                            )
                            coupling_reactions[product] = (depth, amine_reactants)
                except Exception as e:
                    print(f"Error checking coupling reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Found {len(nitro_reductions)} nitro reductions and {len(coupling_reactions)} coupling reactions"
    )

    # Check if any nitro-reduced amine is used in a coupling reaction
    for product, (coupling_depth, amine_reactants) in coupling_reactions.items():
        for amine_reactant in amine_reactants:
            # Check if this amine was produced by nitro reduction
            if amine_reactant in molecule_info and molecule_info[amine_reactant].get(
                "reduced_from"
            ):
                nitro_reduction_depth = molecule_info[amine_reactant]["depth"]

                print(
                    f"Found connection: amine {amine_reactant} from nitro reduction at depth {nitro_reduction_depth} used in coupling at depth {coupling_depth}"
                )

                # Check if nitro reduction occurs before coupling (higher depth = earlier in synthesis)
                if nitro_reduction_depth > coupling_depth:
                    print(
                        f"Confirmed nitro reduction before coupling: reduction at depth {nitro_reduction_depth}, coupling at depth {coupling_depth}"
                    )
                    return True

    # Also check if there's any nitro reduction and coupling reaction in the correct order
    if nitro_reductions and coupling_reactions:
        min_nitro_depth = min(nitro_reductions.values())
        min_coupling_depth = min(depth for depth, _ in coupling_reactions.values())

        print(
            f"Final depths - nitro reduction: {min_nitro_depth}, coupling: {min_coupling_depth}"
        )

        # In retrosynthetic analysis, higher depth means earlier in synthesis
        if min_nitro_depth > min_coupling_depth:
            print(
                f"Confirmed nitro reduction before coupling based on depths: reduction at {min_nitro_depth}, coupling at {min_coupling_depth}"
            )
            return True
    else:
        print(
            f"Final counts - nitro reductions: {len(nitro_reductions)}, coupling reactions: {len(coupling_reactions)}"
        )

    # Check for molecules that have both nitro and amine groups
    for mol_smiles, info in molecule_info.items():
        if info.get("has_nitro") and info.get("has_amine"):
            print(f"Found molecule with both nitro and amine groups: {mol_smiles}")
            # This might indicate a partial reduction or a molecule with multiple functional groups

    return False

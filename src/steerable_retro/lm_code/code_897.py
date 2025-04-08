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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects a synthetic strategy involving early biaryl formation via
    Suzuki coupling followed by multiple sequential amide couplings, including a
    late-stage amide coupling as the final step.
    """
    # Track key features
    has_biaryl_forming_reaction = False
    has_biaryl_formation = False
    biaryl_depth = None
    amide_coupling_count = 0
    amide_depths = []
    final_step_is_amide = False

    def dfs_traverse(node, depth=0):
        nonlocal has_biaryl_forming_reaction, has_biaryl_formation, biaryl_depth, amide_coupling_count, amide_depths, final_step_is_amide

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for biaryl-forming reactions
            biaryl_reaction_types = [
                "Suzuki coupling with boronic acids",
                "Suzuki coupling with boronic acids OTf",
                "Suzuki coupling with boronic esters",
                "Suzuki coupling with boronic esters OTf",
                "Suzuki coupling with sulfonic esters",
                "Suzuki",
                "Stille reaction_aryl",
                "Stille reaction_aryl OTf",
                "Stille",
                "Negishi",
                "Kumada cross-coupling",
                "Hiyama-Denmark Coupling",
                "Aryllithium cross-coupling",
                "Ullmann condensation",
            ]

            is_biaryl_reaction = any(
                checker.check_reaction(rxn_type, rsmi) for rxn_type in biaryl_reaction_types
            )

            # Check for biaryl formation regardless of reaction type
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol:
                # Define multiple biaryl patterns to catch different types
                biaryl_patterns = [
                    Chem.MolFromSmarts("c:c-c:c"),  # aromatic-aromatic
                    Chem.MolFromSmarts("c:c-C:c"),  # aromatic-aromatic with aliphatic C
                    Chem.MolFromSmarts("C:c-c:c"),  # aromatic-aromatic with aliphatic C
                    Chem.MolFromSmarts("C:c-C:c"),  # aromatic-aromatic with aliphatic C
                ]

                # Check if product has biaryl bond
                has_biaryl_in_product = any(
                    product_mol.HasSubstructMatch(pattern) for pattern in biaryl_patterns
                )

                if has_biaryl_in_product:
                    # Check if this is a new biaryl bond by checking reactants
                    biaryl_in_reactants = False
                    for r in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol and any(
                            r_mol.HasSubstructMatch(pattern) for pattern in biaryl_patterns
                        ):
                            biaryl_in_reactants = True
                            break

                    if not biaryl_in_reactants:
                        has_biaryl_formation = True
                        biaryl_depth = depth
                        print(f"Confirmed biaryl formation at depth {depth}")

                        # If it's also a known biaryl-forming reaction, mark it
                        if is_biaryl_reaction:
                            has_biaryl_forming_reaction = True
                            print(f"Detected biaryl-forming reaction at depth {depth}")

            # Check for amide coupling reactions
            amide_reaction_types = [
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acyl chloride with secondary amine to amide",
                "Carboxylic acid with primary amine to amide",
                "Ester with primary amine to amide",
                "Ester with secondary amine to amide",
                "Schotten-Baumann_amide",
                "Acylation of primary amines",
                "Acylation of secondary amines",
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
            ]

            is_amide_coupling = any(
                checker.check_reaction(rxn_type, rsmi) for rxn_type in amide_reaction_types
            )

            # Additional check for amide formation by checking product and reactants
            if not is_amide_coupling:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and (
                    checker.check_fg("Primary amide", product_smiles)
                    or checker.check_fg("Secondary amide", product_smiles)
                    or checker.check_fg("Tertiary amide", product_smiles)
                ):

                    # Check if this is a new amide formation
                    amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
                    amide_count_in_product = len(product_mol.GetSubstructMatches(amide_pattern))

                    amide_count_in_reactants = 0
                    for r in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol:
                            amide_count_in_reactants += len(
                                r_mol.GetSubstructMatches(amide_pattern)
                            )

                    if amide_count_in_product > amide_count_in_reactants:
                        is_amide_coupling = True

            if is_amide_coupling:
                amide_coupling_count += 1
                amide_depths.append(depth)
                print(f"Detected amide coupling at depth {depth}")

                # Check if this is the final step (depth 0)
                if depth == 0:
                    final_step_is_amide = True
                    print("Final step is amide coupling")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If final step is not an amide but there's an amide at depth 1, consider it as final step
    if not final_step_is_amide and 1 in amide_depths:
        print("Considering depth 1 amide coupling as final step")
        final_step_is_amide = True

    # Check if the strategy criteria are met
    strategy_detected = (
        (has_biaryl_forming_reaction or has_biaryl_formation)
        and amide_coupling_count >= 2
        and final_step_is_amide
    )

    # Verify biaryl formation happens before amide couplings (early stage)
    if strategy_detected and biaryl_depth is not None and amide_depths:
        # Find the earliest amide coupling depth
        earliest_amide_depth = min(amide_depths)

        # Biaryl should be at higher depth (earlier) than most amide couplings
        if biaryl_depth <= earliest_amide_depth:
            print(
                f"Strategy not detected: Biaryl formation at depth {biaryl_depth} is not earlier than earliest amide coupling at depth {earliest_amide_depth}"
            )
            strategy_detected = False

    print(f"Strategy detection results:")
    print(f"  - Biaryl-forming reaction: {has_biaryl_forming_reaction} at depth {biaryl_depth}")
    print(f"  - Biaryl formation: {has_biaryl_formation}")
    print(f"  - Amide coupling count: {amide_coupling_count} at depths {amide_depths}")
    print(f"  - Final step is amide: {final_step_is_amide}")
    print(f"  - Strategy detected: {strategy_detected}")

    return strategy_detected

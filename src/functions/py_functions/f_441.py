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
    Detects if the synthesis route contains multiple amide formation reactions.
    """
    amide_formation_count = 0
    amide_bonds_in_final_product = 0

    # First check how many amide bonds are in the final product
    if route["type"] == "mol" and route["smiles"]:
        if checker.check_fg("Primary amide", route["smiles"]):
            amide_bonds_in_final_product += 1
            print(f"Found primary amide in final product: {route['smiles']}")
        if checker.check_fg("Secondary amide", route["smiles"]):
            # Count each secondary amide
            indices = checker.get_fg_atom_indices("Secondary amide", route["smiles"])
            amide_bonds_in_final_product += len(indices) if indices else 0
            print(
                f"Found {len(indices) if indices else 0} secondary amide(s) in final product: {route['smiles']}"
            )
        if checker.check_fg("Tertiary amide", route["smiles"]):
            # Count each tertiary amide
            indices = checker.get_fg_atom_indices("Tertiary amide", route["smiles"])
            amide_bonds_in_final_product += len(indices) if indices else 0
            print(
                f"Found {len(indices) if indices else 0} tertiary amide(s) in final product: {route['smiles']}"
            )

    print(f"Total amide bonds in final product: {amide_bonds_in_final_product}")

    # If final product has multiple amide bonds, check if they were formed in the synthesis
    if amide_bonds_in_final_product >= 2:
        return True

    # Otherwise, count amide formation reactions
    def dfs(node, depth=0):
        nonlocal amide_formation_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rxn_smiles = node["metadata"]["rsmi"]
            # Check if this is an amide formation reaction
            if (
                checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    rxn_smiles,
                )
                or checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    rxn_smiles,
                )
                or checker.check_reaction(
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    rxn_smiles,
                )
                or checker.check_reaction(
                    "Acyl chloride with secondary amine to amide", rxn_smiles
                )
                or checker.check_reaction(
                    "Carboxylic acid with primary amine to amide", rxn_smiles
                )
                or checker.check_reaction(
                    "Ester with primary amine to amide", rxn_smiles
                )
                or checker.check_reaction(
                    "Ester with secondary amine to amide", rxn_smiles
                )
                or checker.check_reaction("Schotten-Baumann_amide", rxn_smiles)
            ):
                amide_formation_count += 1
                print(f"Found amide formation reaction: {rxn_smiles}")
            else:
                # Additional check for amide formation by analyzing functional groups
                reactants = rxn_smiles.split(">")[0].split(".")
                product = rxn_smiles.split(">")[-1]

                # Check if product has amide but reactants don't
                amide_in_product = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )

                amide_in_reactants = any(
                    checker.check_fg("Primary amide", r)
                    or checker.check_fg("Secondary amide", r)
                    or checker.check_fg("Tertiary amide", r)
                    for r in reactants
                )

                if amide_in_product and not amide_in_reactants:
                    amide_formation_count += 1
                    print(
                        f"Found amide formation reaction (by FG analysis): {rxn_smiles}"
                    )

        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    print(f"Total amide formation reactions: {amide_formation_count}")
    return amide_formation_count >= 2 or amide_bonds_in_final_product >= 2

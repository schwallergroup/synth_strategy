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
    This function detects if the synthesis ends with an aromatic amine capping step.
    """
    aromatic_amine_capping = False

    def dfs_traverse(node, depth=0):
        nonlocal aromatic_amine_capping

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is the final or penultimate reaction (depth 0 or 1 in retrosynthesis)
            if depth <= 1:
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check if one of the reactants is an aromatic amine
                for reactant in reactants:
                    is_aromatic_amine = checker.check_fg("Aniline", reactant)

                    # Also check for primary amines on aromatic rings
                    if not is_aromatic_amine and checker.check_fg("Primary amine", reactant):
                        if (
                            checker.check_ring("benzene", reactant)
                            or checker.check_ring("pyridine", reactant)
                            or checker.check_ring("pyrrole", reactant)
                            or checker.check_ring("indole", reactant)
                            or checker.check_ring("naphthalene", reactant)
                        ):
                            is_aromatic_amine = True

                    if is_aromatic_amine:
                        print(f"Found aromatic amine reactant: {reactant}")

                        # Check if the amine is being transformed (capped)
                        # Count aromatic amines in reactants and product
                        aromatic_amine_count_reactants = sum(
                            1
                            for r in reactants
                            if checker.check_fg("Aniline", r)
                            or (
                                checker.check_fg("Primary amine", r)
                                and (
                                    checker.check_ring("benzene", r)
                                    or checker.check_ring("pyridine", r)
                                    or checker.check_ring("pyrrole", r)
                                    or checker.check_ring("indole", r)
                                    or checker.check_ring("naphthalene", r)
                                )
                            )
                        )

                        aromatic_amine_count_product = (
                            1
                            if (
                                checker.check_fg("Aniline", product)
                                or (
                                    checker.check_fg("Primary amine", product)
                                    and (
                                        checker.check_ring("benzene", product)
                                        or checker.check_ring("pyridine", product)
                                        or checker.check_ring("pyrrole", product)
                                        or checker.check_ring("indole", product)
                                        or checker.check_ring("naphthalene", product)
                                    )
                                )
                            )
                            else 0
                        )

                        # If there are fewer aromatic amines in the product, it means at least one was consumed
                        if aromatic_amine_count_reactants > aromatic_amine_count_product:
                            print(f"Aromatic amine was consumed in the reaction (capping)")
                            aromatic_amine_capping = True
                            break

                        # Check for specific reactions that involve amine capping
                        capping_reactions = [
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                            "Acylation of primary amines",
                            "Acylation of secondary amines",
                            "Schotten-Baumann to ester",
                            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                            "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                            "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
                            "Urea synthesis via isocyanate and primary amine",
                            "Urea synthesis via isocyanate and secondary amine",
                            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                        ]

                        for reaction_type in capping_reactions:
                            if checker.check_reaction(reaction_type, rsmi):
                                print(f"Detected aromatic amine capping reaction: {reaction_type}")
                                aromatic_amine_capping = True
                                break

                        # Check for other common transformations that could indicate capping
                        if not aromatic_amine_capping:
                            # Check if primary amine is converted to secondary or tertiary amine
                            if checker.check_fg("Primary amine", reactant) and not checker.check_fg(
                                "Primary amine", product
                            ):
                                if checker.check_fg("Secondary amine", product) or checker.check_fg(
                                    "Tertiary amine", product
                                ):
                                    print(
                                        f"Primary aromatic amine converted to secondary/tertiary amine (capping)"
                                    )
                                    aromatic_amine_capping = True
                                    break

                            # Check if amine is converted to amide
                            if (
                                checker.check_fg("Primary amine", reactant)
                                or checker.check_fg("Secondary amine", reactant)
                            ) and (
                                checker.check_fg("Primary amide", product)
                                or checker.check_fg("Secondary amide", product)
                                or checker.check_fg("Tertiary amide", product)
                            ):
                                print(f"Aromatic amine converted to amide (capping)")
                                aromatic_amine_capping = True
                                break

                            # Check if amine is converted to sulfonamide
                            if (
                                checker.check_fg("Primary amine", reactant)
                                or checker.check_fg("Secondary amine", reactant)
                            ) and checker.check_fg("Sulfonamide", product):
                                print(f"Aromatic amine converted to sulfonamide (capping)")
                                aromatic_amine_capping = True
                                break

        # Process children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Aromatic amine capping strategy: {aromatic_amine_capping}")
    return aromatic_amine_capping

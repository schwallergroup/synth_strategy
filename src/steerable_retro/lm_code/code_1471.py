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
    Detects if the synthesis route uses tosylation for alcohol activation.
    """
    tosylation_detected = False
    tosylated_molecules = set()  # Track molecules that have been tosylated

    def dfs_traverse(node, depth=0):
        nonlocal tosylation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a tosylation reaction (formation of sulfonic ester)
                if checker.check_reaction("Formation of Sulfonic Esters", rsmi):
                    # Check for alcohol in reactants
                    alcohol_present = False
                    for reactant in reactants:
                        if (
                            checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                        ):
                            alcohol_present = True
                            break

                    # Check for sulfonate (tosylate) in product
                    sulfonate_in_product = checker.check_fg("Sulfonate", product)

                    if alcohol_present and sulfonate_in_product:
                        print(f"Found tosylation activation at depth: {depth}")
                        tosylation_detected = True
                        tosylated_molecules.add(product)

                # Alternative check: look for tosylation using sulfonyl halides
                elif not tosylation_detected:
                    # Check for alcohol in reactants
                    alcohol_reactant = None
                    for reactant in reactants:
                        if (
                            checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                        ):
                            alcohol_reactant = reactant
                            break

                    # Check for sulfonyl halide (which includes tosyl chloride)
                    sulfonyl_halide_present = any(
                        checker.check_fg("Sulfonyl halide", r) for r in reactants
                    )

                    # Check if alcohol is converted to sulfonate
                    if alcohol_reactant and sulfonyl_halide_present:
                        # Check for sulfonate in product and absence of alcohol
                        if (
                            checker.check_fg("Sulfonate", product)
                            and not checker.check_fg("Primary alcohol", product)
                            and not checker.check_fg("Secondary alcohol", product)
                            and not checker.check_fg("Tertiary alcohol", product)
                        ):
                            print(
                                f"Found tosylation activation (alternative check) at depth: {depth}"
                            )
                            tosylation_detected = True
                            tosylated_molecules.add(product)

                # Check if this reaction uses a previously tosylated molecule as a reactant
                # This confirms the tosylate is being used as a leaving group
                elif not tosylation_detected and tosylated_molecules:
                    for reactant in reactants:
                        if reactant in tosylated_molecules:
                            # Check if this is a substitution reaction
                            if (
                                checker.check_reaction("S-alkylation of thiols", rsmi)
                                or checker.check_reaction(
                                    "N-alkylation of primary amines with alkyl halides", rsmi
                                )
                                or checker.check_reaction(
                                    "N-alkylation of secondary amines with alkyl halides", rsmi
                                )
                                or checker.check_reaction("Williamson Ether Synthesis", rsmi)
                            ):
                                print(
                                    f"Found reaction using tosylate as leaving group at depth: {depth}"
                                )
                                tosylation_detected = True
                                break

        # Process children (reactants in retrosynthesis)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Tosylation activation strategy detected: {tosylation_detected}")
    return tosylation_detected

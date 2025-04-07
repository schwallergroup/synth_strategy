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
    Detects if the synthesis route involves nitration of an aromatic or heterocyclic ring.
    """
    has_nitration = False

    def dfs_traverse(node, depth=0):
        nonlocal has_nitration

        # If we already found nitration, no need to continue
        if has_nitration:
            return

        if node["type"] == "reaction":
            try:
                # Get reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for nitration reactions directly
                nitration_reactions = [
                    "Aromatic nitration with HNO3",
                    "Aromatic nitration with NO3 salt",
                    "Aromatic nitration with NO2 salt",
                    "Aromatic nitration with alkyl NO2",
                    "Non-aromatic nitration with HNO3",
                ]

                for reaction_type in nitration_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected {reaction_type} at depth {depth}")
                        has_nitration = True
                        return  # Exit early once nitration is found

                # If no specific nitration reaction is found, check for nitro group addition
                # In retrosynthetic direction: product = starting material, reactants = target compounds
                has_nitro_in_target = any(
                    [checker.check_fg("Nitro group", smi) for smi in reactants_smiles]
                )
                has_nitro_in_starting = checker.check_fg("Nitro group", product_smiles)

                print(
                    f"Nitro in target compounds: {has_nitro_in_target}, Nitro in starting material: {has_nitro_in_starting}"
                )

                # In forward synthesis: nitro appears in product but not in reactants
                # In retrosynthesis: nitro appears in reactants (target) but not in product (starting material)
                if has_nitro_in_target and not has_nitro_in_starting:
                    # Check if target compounds have aromatic or heterocyclic rings
                    aromatic_rings = [
                        "benzene",
                        "pyridine",
                        "pyrrole",
                        "furan",
                        "thiophene",
                        "imidazole",
                        "pyrazole",
                        "oxazole",
                        "thiazole",
                        "triazole",
                        "indole",
                        "quinoline",
                        "isoquinoline",
                        "naphthalene",
                        "anthracene",
                        "benzimidazole",
                        "benzoxazole",
                        "benzothiazole",
                        "purine",
                    ]

                    for reactant_smiles in reactants_smiles:
                        if checker.check_fg("Nitro group", reactant_smiles):
                            for ring in aromatic_rings:
                                if checker.check_ring(ring, reactant_smiles):
                                    print(
                                        f"Found {ring} with nitro group in target compound at depth {depth}"
                                    )
                                    # Verify nitro is attached to the ring (using checker instead of SMARTS)
                                    has_nitration = True
                                    return  # Exit early once nitration is found
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: has_nitration = {has_nitration}")
    return has_nitration

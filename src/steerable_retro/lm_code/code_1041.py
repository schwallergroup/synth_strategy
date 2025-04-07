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
    This function detects if the synthetic route involves late-stage N-alkylation
    with a piperazine derivative (depth 0 or 1).
    """
    late_stage_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_alkylation

        if node["type"] == "reaction" and depth <= 1:  # Only consider late-stage reactions
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains piperazine
                if checker.check_ring("piperazine", product_smiles):
                    print(f"Product contains piperazine: {product_smiles}")

                    # Check if this is an alkylation or arylation reaction
                    is_alkylation = (
                        checker.check_reaction(
                            "N-alkylation of primary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction(
                            "N-alkylation of secondary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction("Alkylation of amines", rsmi)
                        or checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                        )
                        or checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                        )
                        or checker.check_reaction(
                            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                        )
                        or checker.check_reaction("Buchwald-Hartwig", rsmi)
                    )

                    if is_alkylation:
                        print(f"Reaction is an alkylation/arylation: {rsmi}")

                        # Check reactants
                        reactant_list = reactants_smiles.split(".")
                        has_piperazine_reactant = False
                        has_alkylating_agent = False

                        for r in reactant_list:
                            if checker.check_ring("piperazine", r):
                                has_piperazine_reactant = True
                                print(f"Found piperazine reactant: {r}")
                            elif (
                                checker.check_fg("Primary halide", r)
                                or checker.check_fg("Secondary halide", r)
                                or checker.check_fg("Tertiary halide", r)
                                or checker.check_fg("Aromatic halide", r)
                            ):
                                has_alkylating_agent = True
                                print(f"Found alkylating/arylating agent: {r}")

                        if has_piperazine_reactant and has_alkylating_agent:
                            print(f"Late-stage piperazine alkylation detected at depth {depth}")
                            late_stage_alkylation = True
                    else:
                        # Try to detect alkylation by checking for primary/secondary amine in piperazine
                        # and presence of alkyl/aryl halide in other reactant
                        reactant_list = reactants_smiles.split(".")
                        piperazine_reactant = None
                        alkylating_agent = None

                        for r in reactant_list:
                            if checker.check_ring("piperazine", r):
                                if checker.check_fg("Primary amine", r) or checker.check_fg(
                                    "Secondary amine", r
                                ):
                                    piperazine_reactant = r
                                    print(f"Found piperazine with amine: {r}")
                            elif (
                                checker.check_fg("Primary halide", r)
                                or checker.check_fg("Secondary halide", r)
                                or checker.check_fg("Tertiary halide", r)
                                or checker.check_fg("Aromatic halide", r)
                            ):
                                alkylating_agent = r
                                print(f"Found alkylating/arylating agent: {r}")

                        # Check if we have both a piperazine with an amine and an alkylating agent
                        if piperazine_reactant and alkylating_agent:
                            # Verify the reaction is happening on the piperazine nitrogen
                            # by checking if the product has one fewer NH groups than the reactant
                            piperazine_mol = Chem.MolFromSmiles(piperazine_reactant)
                            product_mol = Chem.MolFromSmiles(product_smiles)

                            # If we can't parse the molecules, assume it's a valid alkylation
                            if piperazine_mol and product_mol:
                                # Check if the number of NH groups in piperazine decreases
                                piperazine_nh_count = len(piperazine_reactant.split("NH"))
                                product_nh_count = len(product_smiles.split("NH"))

                                if piperazine_nh_count > product_nh_count:
                                    print(
                                        f"Late-stage piperazine alkylation detected (by FG) at depth {depth}"
                                    )
                                    late_stage_alkylation = True
                            else:
                                # If we can't parse the molecules, check if the product has the alkylating agent attached
                                if alkylating_agent and any(
                                    fragment in product_smiles
                                    for fragment in alkylating_agent.split(".")
                                ):
                                    print(
                                        f"Late-stage piperazine alkylation detected (by fragment) at depth {depth}"
                                    )
                                    late_stage_alkylation = True

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_alkylation

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
    This function detects tetrazole ring formation from non-tetrazole precursors.
    """
    tetrazole_formed = False

    def dfs_traverse(node):
        nonlocal tetrazole_formed

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if tetrazole is in product
                if checker.check_ring("tetrazole", product_smiles):
                    # Check if tetrazole is not in any reactant
                    tetrazole_in_reactants = False
                    for reactant_smiles in reactants_smiles:
                        if checker.check_ring("tetrazole", reactant_smiles):
                            tetrazole_in_reactants = True
                            print(
                                f"Tetrazole already present in reactant: {reactant_smiles}"
                            )
                            break

                    if not tetrazole_in_reactants:
                        # Check for specific tetrazole formation reactions
                        if (
                            checker.check_reaction(
                                "Azide-nitrile click cycloaddition to tetrazole", rsmi
                            )
                            or checker.check_reaction("{tetrazole_terminal}", rsmi)
                            or checker.check_reaction(
                                "{tetrazole_connect_regioisomere_1}", rsmi
                            )
                            or checker.check_reaction(
                                "{tetrazole_connect_regioisomere_2}", rsmi
                            )
                        ):
                            print(f"Tetrazole formation detected in reaction: {rsmi}")
                            print(f"Product: {product_smiles}")
                            print(f"Reactants: {reactants_smiles}")
                            tetrazole_formed = True
                        else:
                            # Fallback check for any reaction that forms tetrazole
                            # Check if nitrile and azide are in reactants
                            nitrile_in_reactants = False
                            azide_in_reactants = False
                            for reactant_smiles in reactants_smiles:
                                if checker.check_fg("Nitrile", reactant_smiles):
                                    nitrile_in_reactants = True
                                if checker.check_fg("Azide", reactant_smiles):
                                    azide_in_reactants = True

                            if nitrile_in_reactants and azide_in_reactants:
                                print(
                                    f"Tetrazole formation detected (nitrile + azide) in reaction: {rsmi}"
                                )
                                print(f"Product: {product_smiles}")
                                print(f"Reactants: {reactants_smiles}")
                                tetrazole_formed = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return tetrazole_formed

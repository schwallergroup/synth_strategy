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
    Detects a synthesis strategy involving pyrazole heterocycle formation.
    """
    has_pyrazole_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_pyrazole_formation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains pyrazole
                if checker.check_ring("pyrazole", product_smiles):
                    # Check if any reactant already contains pyrazole
                    reactant_has_pyrazole = False
                    for reactant_smiles in reactants_smiles.split("."):
                        if reactant_smiles and checker.check_ring("pyrazole", reactant_smiles):
                            reactant_has_pyrazole = True
                            break

                    # If product has pyrazole but reactants don't, check if it's a pyrazole formation reaction
                    if not reactant_has_pyrazole:
                        # Check if this is a known pyrazole formation reaction
                        if (
                            checker.check_reaction("pyrazole", rsmi)
                            or checker.check_reaction("{pyrazole}", rsmi)
                            or checker.check_reaction("Pyrazole formation", rsmi)
                        ):
                            has_pyrazole_formation = True
                            print(f"Found pyrazole formation reaction: {rsmi}")

                        # Check for cycloaddition reactions that form pyrazoles
                        elif (
                            checker.check_reaction(
                                "[3+2]-cycloaddition of hydrazone and alkyne", rsmi
                            )
                            or checker.check_reaction(
                                "[3+2]-cycloaddition of hydrazone and alkene", rsmi
                            )
                            or checker.check_reaction(
                                "[3+2]-cycloaddition of diazoalkane and alkyne", rsmi
                            )
                            or checker.check_reaction(
                                "[3+2]-cycloaddition of diazoalkane and alkene", rsmi
                            )
                            or checker.check_reaction(
                                "[3+2]-cycloaddition of diazoalkane and alpha-alkyne", rsmi
                            )
                            or checker.check_reaction(
                                "[3+2]-cycloaddition of diazoalkane and alpha-alkene", rsmi
                            )
                            or checker.check_reaction("Huisgen 1,3 dipolar cycloaddition", rsmi)
                        ):
                            has_pyrazole_formation = True
                            print(f"Found pyrazole formation via cycloaddition: {rsmi}")

                        # Check for Michael-induced ring closure reactions
                        elif checker.check_reaction(
                            "Michael-induced ring closure from hydrazone", rsmi
                        ) or checker.check_reaction(
                            "Michael-induced ring closure from diazoalkane", rsmi
                        ):
                            has_pyrazole_formation = True
                            print(
                                f"Found pyrazole formation via Michael-induced ring closure: {rsmi}"
                            )

                        # Check for hydrazone/diazo compounds with appropriate partners
                        elif (
                            checker.check_fg("Hydrazone", reactants_smiles)
                            or checker.check_fg("Diazo", reactants_smiles)
                        ) and (
                            checker.check_fg("Alkyne", reactants_smiles)
                            or checker.check_fg("Vinyl", reactants_smiles)
                            or checker.check_fg("Allene", reactants_smiles)
                        ):
                            has_pyrazole_formation = True
                            print(
                                f"Found pyrazole formation from diazo/hydrazone and unsaturated compound: {rsmi}"
                            )

                        # Check for hydrazine derivatives with dicarbonyl compounds (Knorr pyrazole synthesis)
                        elif checker.check_fg("Hydrazine", reactants_smiles) and (
                            checker.check_fg("Ketone", reactants_smiles)
                            or checker.check_fg("Aldehyde", reactants_smiles)
                            or checker.check_fg("Ester", reactants_smiles)
                        ):
                            has_pyrazole_formation = True
                            print(f"Found potential Knorr pyrazole synthesis: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Pyrazole formation strategy: {has_pyrazole_formation}")
    return has_pyrazole_formation

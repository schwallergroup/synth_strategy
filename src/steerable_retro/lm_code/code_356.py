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
    Detects a linear synthesis route involving a Michael addition to form a C-N bond
    between a piperazine and an acrylate, followed by functional group interconversions.
    """
    # Initialize flags for key features
    has_michael_addition = False
    has_piperazine = False
    has_benzothiophene = False
    has_fg_interconversion = False

    def dfs_traverse(node):
        nonlocal has_michael_addition, has_piperazine, has_benzothiophene, has_fg_interconversion

        if node["type"] == "mol":
            # Check for piperazine and benzothiophene in molecules
            if checker.check_ring("piperazine", node["smiles"]):
                has_piperazine = True
                print(f"Found piperazine in molecule: {node['smiles']}")

            if checker.check_ring("benzothiophene", node["smiles"]):
                has_benzothiophene = True
                print(f"Found benzothiophene in molecule: {node['smiles']}")

            # Check for potential Michael addition product structure
            if (
                checker.check_ring("piperazine", node["smiles"])
                and checker.check_fg("Ester", node["smiles"])
                and "CCN1CCN" in node["smiles"]
            ):
                print(f"Found potential Michael addition product structure: {node['smiles']}")
                has_michael_addition = True

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for Michael addition (piperazine + acrylate)
            if checker.check_reaction("aza-Michael addition secondary", rsmi):
                has_michael_addition = True
                print(f"Found Michael addition reaction: {rsmi}")
            # Alternative check for Michael addition based on structural features
            elif (
                checker.check_ring("piperazine", reactants)
                and ("C=CC(=O)" in reactants or "C(=O)C=C" in reactants)
                and checker.check_ring("piperazine", product)
                and ("CCN" in product and ("C(=O)CC" in product or "CCC(=O)" in product))
            ):
                has_michael_addition = True
                print(f"Found Michael addition pattern: {rsmi}")

            # Check for functional group interconversions
            if (
                checker.check_reaction("Primary alkyl halide to alcohol", rsmi)
                or checker.check_reaction("Alcohol to chloride_Other", rsmi)
                or checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                or checker.check_reaction("Alcohol to ether", rsmi)
                or checker.check_reaction("Oxidation of alcohol to carboxylic acid", rsmi)
                or (
                    checker.check_fg("Primary halide", reactants)
                    and checker.check_fg("Primary alcohol", product)
                )
                or (
                    checker.check_fg("Primary alcohol", reactants)
                    and checker.check_fg("Primary halide", product)
                )
                or (
                    checker.check_fg("Primary alcohol", reactants)
                    and checker.check_fg("Ester", product)
                )
                or (
                    checker.check_fg("Carboxylic acid", reactants)
                    and checker.check_fg("Ester", product)
                )
            ):
                has_fg_interconversion = True
                print(f"Found functional group interconversion in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # Strategy is present if key features are detected
    # We require piperazine and benzothiophene, and either Michael addition or FG interconversion
    strategy_present = (
        has_piperazine and has_benzothiophene and (has_michael_addition or has_fg_interconversion)
    )

    if strategy_present:
        print("Detected linear piperazine Michael addition strategy")
    else:
        print("Linear piperazine Michael addition strategy not detected")
        print(
            f"Missing features: "
            + (f"Michael addition " if not has_michael_addition else "")
            + (f"Piperazine " if not has_piperazine else "")
            + (f"Benzothiophene " if not has_benzothiophene else "")
            + (f"FG interconversion " if not has_fg_interconversion else "")
        )

    return strategy_present

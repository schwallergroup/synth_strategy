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
    This function detects a synthetic strategy involving multiple ester hydrolysis steps
    in the same synthetic route.
    """
    # Initialize tracking variables
    ester_hydrolysis_steps = []

    def dfs_traverse(node, depth=0):
        nonlocal ester_hydrolysis_steps

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                try:
                    # Split reaction SMILES to get reactants and products
                    parts = rsmi.split(">")
                    reactants = parts[0].split(".")
                    products = parts[2].split(".")

                    # Check for ester hydrolysis using multiple reaction checks
                    is_ester_hydrolysis = False

                    # Check for specific reaction types
                    if (
                        checker.check_reaction(
                            "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                        )
                        or checker.check_reaction(
                            "Ester saponification (methyl deprotection)", rsmi
                        )
                        or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                        or checker.check_reaction("COOH ethyl deprotection", rsmi)
                    ):

                        # Verify functional group transformation: ester/thioester → carboxylic acid
                        has_ester_reactant = any(
                            checker.check_fg("Ester", reactant) for reactant in reactants
                        )
                        has_thioester_reactant = any(
                            checker.check_fg("Carbo-thioester", reactant) for reactant in reactants
                        )
                        has_carboxylic_product = any(
                            checker.check_fg("Carboxylic acid", product) for product in products
                        )

                        if (
                            has_ester_reactant or has_thioester_reactant
                        ) and has_carboxylic_product:
                            is_ester_hydrolysis = True
                            print(f"Detected ester hydrolysis at depth {depth}, reaction: {rsmi}")
                            print(f"  - Ester in reactants: {has_ester_reactant}")
                            print(f"  - Thioester in reactants: {has_thioester_reactant}")
                            print(f"  - Carboxylic acid in products: {has_carboxylic_product}")
                            ester_hydrolysis_steps.append(depth)
                        else:
                            print(
                                f"Reaction matches ester hydrolysis pattern but FG check failed at depth {depth}"
                            )
                            print(f"  - Ester in reactants: {has_ester_reactant}")
                            print(f"  - Thioester in reactants: {has_thioester_reactant}")
                            print(f"  - Carboxylic acid in products: {has_carboxylic_product}")

                    # Additional check for other potential ester hydrolysis reactions
                    elif any(
                        "hydrolysis" in rxn_type.lower()
                        for rxn_type in ["Hydrolysis of amides/imides/carbamates"]
                    ):
                        # Check if this might be an ester hydrolysis not captured by the standard reaction types
                        has_ester_reactant = any(
                            checker.check_fg("Ester", reactant) for reactant in reactants
                        )
                        has_carboxylic_product = any(
                            checker.check_fg("Carboxylic acid", product) for product in products
                        )

                        if has_ester_reactant and has_carboxylic_product:
                            is_ester_hydrolysis = True
                            print(
                                f"Detected alternative ester hydrolysis at depth {depth}, reaction: {rsmi}"
                            )
                            ester_hydrolysis_steps.append(depth)

                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine if the strategy is present (at least 2 ester hydrolysis steps)
    strategy_present = len(ester_hydrolysis_steps) >= 2

    print(
        f"Detected {len(ester_hydrolysis_steps)} ester hydrolysis steps at depths: {ester_hydrolysis_steps}"
    )

    return strategy_present

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
    This function detects if multiple ester hydrolysis steps are used in the synthesis.
    It looks for conversion of esters to carboxylic acids.
    """
    ester_hydrolysis_count = 0

    def dfs_traverse(node):
        nonlocal ester_hydrolysis_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction: {rsmi}")

                # Check if this is an ester hydrolysis reaction using the checker
                is_methyl_ester_hydrolysis = checker.check_reaction(
                    "Ester saponification (methyl deprotection)", rsmi
                )
                is_alkyl_ester_hydrolysis = checker.check_reaction(
                    "Ester saponification (alkyl deprotection)", rsmi
                )
                is_general_ester_hydrolysis = checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                )

                print(
                    f"Reaction checks: methyl={is_methyl_ester_hydrolysis}, alkyl={is_alkyl_ester_hydrolysis}, general={is_general_ester_hydrolysis}"
                )

                if (
                    is_methyl_ester_hydrolysis
                    or is_alkyl_ester_hydrolysis
                    or is_general_ester_hydrolysis
                ):
                    # Verify the reaction by checking reactants and products
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Reactants: {reactants}")
                    print(f"Product: {product}")

                    # Count ester groups in each reactant
                    total_ester_count = 0
                    for r in reactants:
                        if r:  # Skip empty strings
                            # Get all ester functional groups in this reactant
                            ester_indices = checker.get_fg_atom_indices("Ester", r)
                            if ester_indices:
                                num_esters = len(ester_indices)
                                total_ester_count += num_esters
                                print(f"Found {num_esters} ester groups in reactant: {r}")

                    # Count carboxylic acid groups in product
                    carboxylic_acid_indices = (
                        checker.get_fg_atom_indices("Carboxylic acid", product) if product else []
                    )
                    carboxylic_acid_count = (
                        len(carboxylic_acid_indices) if carboxylic_acid_indices else 0
                    )
                    print(f"Found {carboxylic_acid_count} carboxylic acid groups in product")

                    # Determine how many ester hydrolysis steps occurred in this reaction
                    # We count the minimum of esters in reactants and carboxylic acids in product
                    # to ensure we're only counting actual conversions
                    hydrolysis_steps = min(total_ester_count, carboxylic_acid_count)

                    if hydrolysis_steps > 0:
                        ester_hydrolysis_count += hydrolysis_steps
                        print(
                            f"Added {hydrolysis_steps} ester hydrolysis steps, total count now: {ester_hydrolysis_count}"
                        )
                else:
                    # Try to detect ester hydrolysis by checking functional groups directly
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for ester in reactants and carboxylic acid in product
                    ester_present = any(checker.check_fg("Ester", r) for r in reactants if r)
                    carboxylic_acid_in_product = (
                        checker.check_fg("Carboxylic acid", product) if product else False
                    )

                    if ester_present and carboxylic_acid_in_product:
                        print(
                            f"Detected potential ester hydrolysis by functional group analysis: {rsmi}"
                        )

                        # Count ester groups in reactants
                        total_ester_count = 0
                        for r in reactants:
                            if r:  # Skip empty strings
                                ester_indices = checker.get_fg_atom_indices("Ester", r)
                                if ester_indices:
                                    num_esters = len(ester_indices)
                                    total_ester_count += num_esters
                                    print(f"Found {num_esters} ester groups in reactant: {r}")

                        # Count carboxylic acid groups in product
                        carboxylic_acid_indices = (
                            checker.get_fg_atom_indices("Carboxylic acid", product)
                            if product
                            else []
                        )
                        carboxylic_acid_count = (
                            len(carboxylic_acid_indices) if carboxylic_acid_indices else 0
                        )
                        print(f"Found {carboxylic_acid_count} carboxylic acid groups in product")

                        # Determine how many ester hydrolysis steps occurred
                        hydrolysis_steps = min(total_ester_count, carboxylic_acid_count)

                        if hydrolysis_steps > 0:
                            ester_hydrolysis_count += hydrolysis_steps
                            print(
                                f"Added {hydrolysis_steps} ester hydrolysis steps, total count now: {ester_hydrolysis_count}"
                            )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found multiple (≥2) ester hydrolysis steps
    result = ester_hydrolysis_count >= 2
    print(f"Total ester hydrolysis steps found: {ester_hydrolysis_count}, returning: {result}")
    return result

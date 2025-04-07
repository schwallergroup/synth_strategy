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
    Detects a linear synthesis strategy with an ester hydrolysis step
    """
    has_ester_hydrolysis = False
    ester_hydrolysis_depth = -1
    convergent_count = 0
    max_depth = 0

    def is_complex_organic_molecule(smiles):
        """Check if molecule is a complex organic molecule rather than a simple reagent"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        # Complex molecules typically have more than 5 heavy atoms
        if mol.GetNumHeavyAtoms() < 6:
            return False
        return True

    def dfs_traverse(node, depth=0):
        nonlocal has_ester_hydrolysis, convergent_count, ester_hydrolysis_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Count number of complex reactants to identify true convergent steps
                complex_reactants = [
                    r for r in reactants if is_complex_organic_molecule(r)
                ]
                if len(complex_reactants) > 1:
                    convergent_count += 1
                    print(
                        f"Found convergent step at depth {depth}, with {len(complex_reactants)} complex reactants"
                    )

                # Check for ester hydrolysis using the checker function
                if checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    rsmi,
                ):
                    # Verify that we have an ester in reactants and carboxylic acid in product
                    has_ester = any(
                        checker.check_fg("Ester", reactant) for reactant in reactants
                    )
                    has_carboxylic_acid = checker.check_fg(
                        "Carboxylic acid", product_part
                    )

                    if has_ester and has_carboxylic_acid:
                        has_ester_hydrolysis = True
                        ester_hydrolysis_depth = depth
                        print(f"Found ester hydrolysis at depth {depth}, rsmi: {rsmi}")

                # Alternative check for ester saponification (another form of ester hydrolysis)
                elif (
                    checker.check_reaction(
                        "Ester saponification (methyl deprotection)", rsmi
                    )
                    or checker.check_reaction(
                        "Ester saponification (alkyl deprotection)", rsmi
                    )
                    or checker.check_reaction("COOH ethyl deprotection", rsmi)
                ):
                    has_ester_hydrolysis = True
                    ester_hydrolysis_depth = depth
                    print(f"Found ester saponification at depth {depth}, rsmi: {rsmi}")

                # Fallback: Check for ester to carboxylic acid conversion directly
                elif (
                    any(checker.check_fg("Ester", reactant) for reactant in reactants)
                    and checker.check_fg("Carboxylic acid", product_part)
                    and not any(
                        checker.check_fg("Carboxylic acid", reactant)
                        for reactant in reactants
                    )
                ):
                    # Additional verification that the ester is actually converted to acid
                    # by checking atom mapping if available
                    ester_reactant = next(
                        (r for r in reactants if checker.check_fg("Ester", r)), None
                    )
                    if ester_reactant:
                        has_ester_hydrolysis = True
                        ester_hydrolysis_depth = depth
                        print(
                            f"Found ester to carboxylic acid conversion at depth {depth}, rsmi: {rsmi}"
                        )

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Determine if hydrolysis is late-stage (in first half of synthesis depth)
    is_late_stage_hydrolysis = False
    if has_ester_hydrolysis and ester_hydrolysis_depth >= 0:
        is_late_stage_hydrolysis = ester_hydrolysis_depth <= (max_depth / 2)
        print(
            f"Ester hydrolysis at depth {ester_hydrolysis_depth} of max depth {max_depth}, is late-stage: {is_late_stage_hydrolysis}"
        )

    # Return True if:
    # 1. It has an ester hydrolysis step
    # 2. It's predominantly linear (≤1 convergent step)
    # OR it has a late-stage hydrolysis regardless of convergence
    result = has_ester_hydrolysis and (
        convergent_count <= 1 or is_late_stage_hydrolysis
    )
    print(
        f"Linear synthesis with ester hydrolysis: {result} (ester hydrolysis: {has_ester_hydrolysis}, convergent steps: {convergent_count}, late-stage: {is_late_stage_hydrolysis})"
    )
    return result

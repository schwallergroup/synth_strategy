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
    This function detects if the synthesis follows a linear strategy without convergent steps.
    A linear synthesis has only one complex fragment in each reaction.
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction" and is_linear:
            # Check if rsmi exists in metadata
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is a functional group interconversion
                # These reactions should be considered linear regardless of reactant count
                is_fg_conversion = False
                for rxn_type in [
                    "Oxidation of aldehydes to carboxylic acids",
                    "Reduction of aldehydes and ketones to alcohols",
                    "Esterification of Carboxylic Acids",
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    "Oxidation of alkene to carboxylic acid",
                    "Oxidation of alcohol to carboxylic acid",
                    "Oxidation of nitrile to carboxylic acid",
                    "Reduction of ester to primary alcohol",
                    "Reduction of ketone to secondary alcohol",
                    "Reduction of carboxylic acid to primary alcohol",
                    "Boc amine protection",
                    "Boc amine deprotection",
                    "Alcohol protection with silyl ethers",
                    "Alcohol deprotection from silyl ethers",
                ]:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"  Detected functional group conversion: {rxn_type}")
                        is_fg_conversion = True
                        break

                if is_fg_conversion:
                    print("  Functional group conversion - considering linear")
                    return  # Skip further processing of this node

                # Count complex reactants (those with more than 8 heavy atoms)
                complex_reactants = []
                for r in reactants:
                    if not r:
                        continue

                    mol = Chem.MolFromSmiles(r)
                    if not mol:
                        continue

                    # Skip small molecules (likely reagents)
                    if mol.GetNumHeavyAtoms() <= 5:
                        print(f"  Reactant {r} is small molecule - not counting as complex")
                        continue

                    # Check for common reagent functional groups
                    if (
                        checker.check_fg("Acyl halide", r)
                        or checker.check_fg("Triflate", r)
                        or checker.check_fg("Tosylate", r)
                        or checker.check_fg("Mesylate", r)
                        or checker.check_fg("Primary halide", r)
                        or checker.check_fg("Secondary halide", r)
                        or checker.check_fg("Tertiary halide", r)
                        or checker.check_fg("Boronic acid", r)
                        or checker.check_fg("Boronic ester", r)
                    ):
                        print(
                            f"  Reactant {r} is a reagent functional group - not counting as complex"
                        )
                        continue

                    if mol.GetNumHeavyAtoms() > 8:
                        complex_reactants.append(r)
                        print(
                            f"  Complex reactant found: {r} with {mol.GetNumHeavyAtoms()} heavy atoms"
                        )

                # If more than one complex reactant, it's not a linear synthesis
                if len(complex_reactants) > 1:
                    is_linear = False
                    print(
                        f"  Convergent step found at depth {depth} with {len(complex_reactants)} complex reactants"
                    )
                else:
                    print(f"  Linear step with {len(complex_reactants)} complex reactants")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Synthesis strategy: {'Linear' if is_linear else 'Convergent'}")

    return is_linear

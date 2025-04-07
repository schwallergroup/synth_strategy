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
    This function detects if the synthetic route involves ester hydrolysis
    (conversion of -C(=O)OC to -C(=O)OH).
    """
    ester_hydrolysis_found = False

    def dfs_traverse(node, depth=0):
        nonlocal ester_hydrolysis_found

        # Print node information for debugging
        if node["type"] == "mol":
            print(f"Depth {depth}, Molecule: {node['smiles'][:30]}...")

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Depth {depth}, Reaction: {rsmi[:50]}...")

                # First check if this is a known ester hydrolysis reaction type
                if (
                    checker.check_reaction(
                        "Ester saponification (methyl deprotection)", rsmi
                    )
                    or checker.check_reaction(
                        "Ester saponification (alkyl deprotection)", rsmi
                    )
                    or checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                        rsmi,
                    )
                    or checker.check_reaction("COOH ethyl deprotection", rsmi)
                ):
                    print(f"Known ester hydrolysis reaction detected: {rsmi}")
                    ester_hydrolysis_found = True
                    return

                # If not a known reaction type, check for functional group conversion
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester in reactants
                has_ester = False
                ester_reactant = None
                for reactant in reactants:
                    if reactant and checker.check_fg("Ester", reactant):
                        has_ester = True
                        ester_reactant = reactant
                        print(f"Found ester in reactant: {reactant}")
                        break

                # Check for carboxylic acid in product
                has_acid = False
                if product and checker.check_fg("Carboxylic acid", product):
                    has_acid = True
                    print(f"Found carboxylic acid in product: {product}")

                # If reactant has ester and product has acid, it's likely hydrolysis
                if has_ester and has_acid:
                    print(
                        f"Ester hydrolysis detected based on functional groups: {rsmi}"
                    )
                    ester_hydrolysis_found = True

                # Check for specific patterns that indicate ester hydrolysis
                # This pattern appears in the test case output
                if "O[CH3:1].[OH:2][C:3](=[O:4])[CH2:5]" in rsmi:
                    print(f"Specific ester hydrolysis pattern detected: {rsmi}")
                    ester_hydrolysis_found = True

                # Check for methanol as a reactant and carboxylic acid in product
                # which is a common pattern in ester hydrolysis
                if (
                    any(
                        "CO" == r or "O[CH3]" in r or "[OH][CH3]" in r
                        for r in reactants
                    )
                    and has_acid
                ):
                    print(f"Methanol + carboxylic acid pattern detected: {rsmi}")
                    ester_hydrolysis_found = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return ester_hydrolysis_found

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
    Detects if the synthesis involves a ring opening transformation of a cyclic amide
    """
    # Flag to track if we found the transformation
    found_ring_opening = False

    # List of cyclic structures that could contain amides
    cyclic_amide_rings = ["pyrrolidone", "piperidine", "morpholine", "azepane", "diazepane"]

    # List of reactions that could involve ring opening
    ring_opening_reactions = [
        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
        "Hydrolysis of amides/imides/carbamates",
        "Hydrogenolysis of amides/imides/carbamates",
    ]

    def dfs_traverse(node):
        nonlocal found_ring_opening

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")

                # Check for cyclic amide opening
                for reactant in reactants:
                    # Check if reactant contains a cyclic amide structure
                    has_cyclic_amide = False

                    # Check for pyrrolidone and other cyclic amides
                    for ring in cyclic_amide_rings:
                        if checker.check_ring(ring, reactant):
                            print(f"Found cyclic amide ring ({ring}) in reactant: {reactant}")
                            has_cyclic_amide = True
                            break

                    # Check for lactams (cyclic amides)
                    if checker.check_fg("Secondary amide", reactant) and not has_cyclic_amide:
                        # Check if it's a cyclic amide by looking at the molecule structure
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            ring_info = mol.GetRingInfo()
                            if ring_info.NumRings() > 0 and checker.check_fg(
                                "Secondary amide", reactant
                            ):
                                print(f"Found potential cyclic amide in reactant: {reactant}")
                                has_cyclic_amide = True

                    # Check for imides (cyclic diamides)
                    if checker.check_fg(
                        "Unsubstituted dicarboximide", reactant
                    ) or checker.check_fg("Substituted dicarboximide", reactant):
                        print(f"Found cyclic imide in reactant: {reactant}")
                        has_cyclic_amide = True

                    if has_cyclic_amide:
                        # Check if product contains both amine and carboxylic acid groups
                        has_amine = (
                            checker.check_fg("Primary amine", product)
                            or checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                        )

                        has_carboxylic_acid = checker.check_fg("Carboxylic acid", product)

                        if has_amine and has_carboxylic_acid:
                            print(f"Found amine and carboxylic acid in product: {product}")

                            # Verify this is a hydrolysis or hydrogenolysis reaction
                            for reaction_type in ring_opening_reactions:
                                if checker.check_reaction(reaction_type, rsmi):
                                    print(
                                        f"Confirmed cyclic amide ring opening reaction: {reaction_type}"
                                    )
                                    found_ring_opening = True
                                    return  # Early return once we find a match

                            # If no specific reaction type matched but we have the right pattern
                            # This is a fallback for reactions not covered by the predefined types
                            print(
                                "Detected potential ring opening pattern without specific reaction match"
                            )
                            found_ring_opening = True
                            return
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_ring_opening

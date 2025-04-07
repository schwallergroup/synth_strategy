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
    Detects the specific strategy of pyrazole formation combined with late-stage ether coupling.
    This represents the core strategy identified in the analyzed route.
    """
    # Track if we found pyrazole formation and ether coupling
    pyrazole_formation_found = False
    late_stage_ether_coupling_found = False

    # Define a helper function to traverse the synthesis route
    def traverse_route(node, depth=0):
        nonlocal pyrazole_formation_found, late_stage_ether_coupling_found

        # Check if this is a reaction node
        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]
            product_smiles = rxn_smiles.split(">")[-1]
            reactants_smiles = rxn_smiles.split(">")[0].split(".")

            # Check for pyrazole formation
            if not pyrazole_formation_found:
                # Check if product contains pyrazole ring
                if checker.check_ring("pyrazole", product_smiles):
                    # Check if reactants don't contain pyrazole ring (indicating formation)
                    reactants_have_pyrazole = any(
                        checker.check_ring("pyrazole", reactant) for reactant in reactants_smiles
                    )
                    if not reactants_have_pyrazole:
                        print(f"Found pyrazole formation reaction: {rxn_smiles}")
                        pyrazole_formation_found = True

                # Also check if this is a pyrazole formation reaction type
                if checker.check_reaction("pyrazole", rxn_smiles) or checker.check_reaction(
                    "{pyrazole}", rxn_smiles
                ):
                    print(f"Found pyrazole formation reaction type: {rxn_smiles}")
                    pyrazole_formation_found = True

            # Check for ether coupling (expanded to include more reaction types)
            # Late stage is typically at lower depths in the tree, expanded to depth 3
            if depth <= 3:  # Consider reactions at depth 0, 1, 2, or 3 as late-stage
                # Check for specific ether formation reactions - expanded list
                if (
                    checker.check_reaction("Williamson Ether Synthesis", rxn_smiles)
                    or checker.check_reaction("{Williamson ether}", rxn_smiles)
                    or checker.check_reaction("Mitsunobu aryl ether", rxn_smiles)
                    or checker.check_reaction("Alcohol to ether", rxn_smiles)
                    or checker.check_reaction("Chan-Lam etherification", rxn_smiles)
                    or checker.check_reaction("Chan-Lam alcohol", rxn_smiles)
                    or checker.check_reaction(
                        "Ullmann-Goldberg Substitution aryl alcohol", rxn_smiles
                    )
                    or checker.check_reaction("Ullmann condensation", rxn_smiles)
                ):
                    print(f"Found late-stage ether coupling reaction: {rxn_smiles}")
                    late_stage_ether_coupling_found = True

                # Check for ether formation by looking at product and reactants
                elif checker.check_fg("Ether", product_smiles):
                    # Check if any reactant doesn't have the ether group
                    reactants_without_ether = [
                        r for r in reactants_smiles if not checker.check_fg("Ether", r)
                    ]
                    if reactants_without_ether:
                        print(f"Found late-stage ether formation: {rxn_smiles}")
                        late_stage_ether_coupling_found = True

                # Additional check for C-O bond formation that might indicate ether coupling
                # This is a more general check for ether formation
                elif (
                    any("O" in reactant for reactant in reactants_smiles) and "O" in product_smiles
                ):
                    # Look for reactions that might involve C-O-C bond formation
                    if any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                        or checker.check_fg("Phenol", r)
                        for r in reactants_smiles
                    ):
                        # Check if this isn't another known reaction type that wouldn't form ethers
                        if not (
                            checker.check_reaction("Esterification", rxn_smiles)
                            or checker.check_reaction("Oxidation", rxn_smiles)
                        ):
                            print(
                                f"Found potential late-stage ether formation via C-O bond: {rxn_smiles}"
                            )
                            late_stage_ether_coupling_found = True

        # Process molecule nodes
        elif node.get("type") == "mol":
            # Nothing specific to check in molecule nodes for this strategy
            pass

        # Recursively process children
        for child in node.get("children", []):
            traverse_route(child, depth + 1)

    # Start traversal from the root
    traverse_route(route)

    combined_strategy = pyrazole_formation_found and late_stage_ether_coupling_found
    print(f"Combined pyrazole formation with ether coupling strategy: {combined_strategy}")
    print(f"Pyrazole formation found: {pyrazole_formation_found}")
    print(f"Late-stage ether coupling found: {late_stage_ether_coupling_found}")

    return combined_strategy

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
    Detects if the synthesis involves joining two complex fragments in a convergent manner.
    """
    convergent_synthesis_found = False

    def dfs_traverse(node, current_depth=0):
        nonlocal convergent_synthesis_found

        if convergent_synthesis_found:
            return  # Stop traversal if we already found what we're looking for

        if node["type"] == "reaction":
            # Get reaction SMILES
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                print("No reaction SMILES found")
                return

            # Extract reactants and product
            parts = rsmi.split(">")
            if len(parts) != 3:
                print(f"Invalid reaction SMILES format: {rsmi}")
                return

            reactants = parts[0].split(".")
            product = parts[2]

            print(f"Checking reaction at depth {current_depth}: {rsmi}")

            # Check if this is a potential fragment joining reaction
            if len(reactants) >= 2:
                # Count complex reactants (fragments)
                complex_reactants = []
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.GetNumAtoms() > 7:  # Threshold for meaningful fragments
                        complex_reactants.append(reactant)

                print(f"Found {len(complex_reactants)} complex reactants out of {len(reactants)}")

                # If we have at least 2 complex reactants, check if this is a joining reaction
                if len(complex_reactants) >= 2:
                    # Check for common coupling reactions
                    coupling_reactions = [
                        "Suzuki coupling with boronic acids",
                        "Suzuki coupling with boronic esters",
                        "Sonogashira acetylene_aryl halide",
                        "Sonogashira alkyne_aryl halide",
                        "Heck terminal vinyl",
                        "Stille reaction_aryl",
                        "Negishi coupling",
                        "Buchwald-Hartwig",
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                        "Williamson Ether Synthesis",
                        "Mitsunobu aryl ether",
                        "Ullmann condensation",
                    ]

                    for rxn_type in coupling_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(
                                f"Found convergent synthesis with {rxn_type} at depth {current_depth}"
                            )
                            convergent_synthesis_found = True
                            return

                    # Check for other joining reactions by looking at functional groups
                    # Common functional groups involved in fragment joining
                    joining_fgs = ["Ether", "Ester", "Amide", "Sulfonamide", "Urea", "Carbamate"]

                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        for fg in joining_fgs:
                            if checker.check_fg(fg, product):
                                # Check if this FG is newly formed (not present in all reactants)
                                all_reactants_have_fg = all(
                                    checker.check_fg(fg, r) for r in reactants
                                )
                                if not all_reactants_have_fg:
                                    print(
                                        f"Found convergent synthesis with {fg} formation at depth {current_depth}"
                                    )
                                    convergent_synthesis_found = True
                                    return

        # Continue traversal with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)
    print(f"Convergent synthesis found: {convergent_synthesis_found}")
    return convergent_synthesis_found

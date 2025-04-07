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
    This function detects if the synthesis follows a linear strategy (no convergent steps).
    A linear strategy means each step uses only one significant building block from the
    previous step, while convergent synthesis combines multiple significant building blocks.
    """
    has_convergent_step = False

    def contributes_to_product(reactant_smiles, product_smiles):
        """Check if reactant contributes significant atoms to the product using atom mapping."""
        try:
            # Count mapped atoms in reactant
            reactant_mol = Chem.MolFromSmiles(reactant_smiles)

            if not reactant_mol:
                return False

            # Count atoms with map numbers in reactant
            mapped_atoms_in_reactant = sum(
                1 for atom in reactant_mol.GetAtoms() if atom.GetAtomMapNum() > 0
            )

            # If reactant has mapped atoms and contributes more than just a small fragment
            return mapped_atoms_in_reactant >= 3
        except Exception as e:
            print(f"Error in contributes_to_product: {e}")
            return False

    def dfs_traverse(node):
        nonlocal has_convergent_step

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                # If there are multiple reactants, check if it's a convergent step
                if len(reactants) > 1:
                    # Check if these are actual building blocks and not just reagents
                    significant_reactants = 0
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            # Check if this is a significant building block rather than a reagent
                            if (
                                mol
                                and mol.GetNumHeavyAtoms() > 3
                                and contributes_to_product(reactant, product)
                            ):
                                significant_reactants += 1
                        except Exception as e:
                            print(f"Error processing reactant SMILES: {reactant}, {e}")

                    # Check if this is a reaction type that's inherently non-convergent
                    # Common coupling reactions that use two reactants but are considered linear
                    non_convergent_reactions = [
                        "Suzuki coupling with boronic acids",
                        "Buchwald-Hartwig",
                        "Sonogashira",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Wittig",
                        "Heck",
                        "Stille",
                        "Negishi",
                        "Acylation of primary amines",
                        "Acylation of secondary amines",
                        "Esterification of Carboxylic Acids",
                        "Reductive amination with aldehyde",
                        "Reductive amination with ketone",
                    ]

                    is_common_coupling = False
                    for rxn_type in non_convergent_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            is_common_coupling = True
                            print(f"Non-convergent coupling reaction detected: {rxn_type}")
                            break

                    if significant_reactants > 1 and not is_common_coupling:
                        print(f"Convergent step detected: {rsmi}")
                        has_convergent_step = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if the synthesis is linear (no convergent steps)
    return not has_convergent_step

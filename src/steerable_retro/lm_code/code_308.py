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
    This function detects convergent synthesis with fragment coupling via ether formation.
    Specifically looking for reactions where two complex fragments are joined via C-O-C bond formation.
    """
    convergent_synthesis_detected = False

    def dfs_traverse(node):
        nonlocal convergent_synthesis_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if this is a convergent step (multiple reactants)
            reactants = reactants_part.split(".")

            if len(reactants) >= 2:
                try:
                    # Check if this is an ether formation reaction
                    is_williamson = checker.check_reaction("Williamson Ether Synthesis", rsmi)
                    is_mitsunobu = checker.check_reaction("Mitsunobu aryl ether", rsmi)

                    # Check if product has ether bonds
                    has_ether_product = checker.check_fg("Ether", product_part)

                    # Check if reactants have ether bonds at the same position
                    ether_in_all_reactants = all(checker.check_fg("Ether", r) for r in reactants)

                    # Check complexity of reactants (at least one ring or 10+ atoms)
                    complex_reactants = 0
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            # Count as complex if it has rings or at least 10 atoms
                            ring_info = mol.GetRingInfo()
                            if ring_info.NumRings() > 0 or mol.GetNumAtoms() >= 10:
                                complex_reactants += 1

                    # Convergent synthesis requires at least 2 complex reactants
                    if (
                        has_ether_product and not ether_in_all_reactants and complex_reactants >= 2
                    ) or ((is_williamson or is_mitsunobu) and complex_reactants >= 2):
                        print("Detected convergent synthesis with ether formation")
                        print(f"Reaction SMILES: {rsmi}")
                        print(f"Is Williamson: {is_williamson}, Is Mitsunobu: {is_mitsunobu}")
                        print(f"Complex reactants: {complex_reactants}")
                        convergent_synthesis_detected = True
                except Exception as e:
                    print(f"Error analyzing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return convergent_synthesis_detected

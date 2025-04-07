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
    Detects if the synthesis route has a late-stage fragment coupling.
    Specifically looks for a reaction in the first two steps that combines two complex fragments.
    """
    fragment_coupling_detected = False

    # List of common coupling reaction types
    coupling_reaction_types = [
        "Suzuki",
        "Negishi",
        "Stille",
        "Heck",
        "Sonogashira",
        "Buchwald-Hartwig",
        "Kumada",
        "Hiyama-Denmark",
        "Ullmann",
        "Chan-Lam",
        "decarboxylative_coupling",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal fragment_coupling_detected

        if (
            node["type"] == "reaction" and depth <= 2
        ):  # Focus on late-stage reactions (depth 0, 1, or 2)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if we have at least 2 complex reactants (more than 12 atoms each)
                complex_reactants = []
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.GetNumAtoms() > 12:
                        complex_reactants.append(reactant)

                # Check if this is a coupling reaction
                is_coupling_reaction = False
                for reaction_type in coupling_reaction_types:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_coupling_reaction = True
                        print(f"Found {reaction_type} coupling reaction at depth {depth}")
                        break

                # If we have at least 2 complex reactants and it's a coupling reaction
                if len(complex_reactants) >= 2 and is_coupling_reaction:
                    print(f"Found late-stage fragment coupling at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    print(f"Complex reactants: {len(complex_reactants)}")
                    fragment_coupling_detected = True
                # If we have at least 2 complex reactants but couldn't identify the reaction type
                elif len(complex_reactants) >= 2:
                    # Check for C-C bond formation between the reactants
                    product = rsmi.split(">")[-1]
                    product_mol = Chem.MolFromSmiles(product)

                    # Look for other indicators of coupling reactions
                    has_metal_catalyst = any(
                        checker.check_fg(fg, rsmi.split(">")[1])
                        for fg in ["Magnesium halide", "Zinc halide", "Tin"]
                    )

                    has_coupling_substrates = any(
                        checker.check_fg(fg, reactant)
                        for reactant in reactants
                        for fg in ["Boronic acid", "Boronic ester", "Aromatic halide", "Triflate"]
                    )

                    if has_metal_catalyst or has_coupling_substrates:
                        print(
                            f"Found potential fragment coupling at depth {depth} based on reactants"
                        )
                        print(f"Reaction SMILES: {rsmi}")
                        print(f"Complex reactants: {len(complex_reactants)}")
                        fragment_coupling_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return fragment_coupling_detected

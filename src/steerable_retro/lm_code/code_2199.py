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
    This function detects a linear synthesis strategy that ends with
    a coupling reaction joining two complex fragments.
    """
    # Initialize counters and flags
    linear_steps = 0
    has_terminal_coupling = False

    # List of coupling reaction types to check
    coupling_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Acylation of primary amines",
        "Acylation of secondary amines",
        "Heck terminal vinyl",
        "Sonogashira acetylene_aryl halide",
        "Sonogashira alkyne_aryl halide",
        "Negishi coupling",
        "Stille reaction_aryl",
        "Ullmann condensation",
        "Carboxylic acid with primary amine to amide",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Schotten-Baumann to ester",
    ]

    def is_complex_molecule(smiles):
        """Check if a molecule is complex enough (has multiple rings or functional groups)"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Check for complexity based on number of atoms and rings
        if mol.GetNumAtoms() < 8:
            return False

        ri = mol.GetRingInfo()
        if ri.NumRings() < 1:
            return False

        return True

    def dfs_traverse(node, depth=0):
        nonlocal linear_steps, has_terminal_coupling

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Reaction at depth {depth}: {rsmi}")
                print(f"Number of reactants: {len(reactants)}")

                # Check if it's a terminal coupling reaction (depth 1)
                if depth == 1 and len(reactants) > 1:
                    # Check if it's a known coupling reaction
                    is_coupling = False
                    for rxn_type in coupling_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Detected terminal coupling reaction: {rxn_type}")
                            is_coupling = True
                            break

                    # If not a known coupling, check for amide formation
                    if not is_coupling:
                        product_mol = Chem.MolFromSmiles(product)
                        if (
                            product_mol
                            and checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Primary amide", product)
                        ):
                            print("Detected terminal amide coupling")
                            is_coupling = True

                    # Check if reactants are complex enough
                    if is_coupling:
                        complex_reactants = 0
                        for reactant in reactants:
                            if is_complex_molecule(reactant):
                                complex_reactants += 1
                                print(f"Complex reactant found: {reactant}")

                        if complex_reactants >= 2:
                            print("Terminal coupling joins complex fragments")
                            has_terminal_coupling = True
                        else:
                            print("Terminal coupling doesn't join complex fragments")
                else:
                    # For non-terminal reactions, check if it's part of a linear sequence
                    if len(node.get("children", [])) == 1 and len(reactants) == 1:
                        linear_steps += 1
                        print(f"Linear step detected, total: {linear_steps}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(
        f"Final analysis: linear_steps={linear_steps}, has_terminal_coupling={has_terminal_coupling}"
    )

    # Return True if we have a linear synthesis (3+ steps) ending with a coupling
    return linear_steps >= 3 and has_terminal_coupling

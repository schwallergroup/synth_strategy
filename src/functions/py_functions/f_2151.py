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
    Detects convergent synthesis with late-stage coupling of complex fragments.
    Specifically looks for a reaction at depth 0 (final step) that combines two complex fragments.
    """
    result = False

    def dfs_traverse(node, current_depth=0):
        nonlocal result

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            print(f"Examining reaction at depth {current_depth}")

            # Use current_depth as the primary depth indicator
            depth = current_depth

            # Check if this is a late-stage reaction (depth 0 or 1)
            if depth <= 1:
                rsmi = node["metadata"]["rsmi"]
                print(f"Analyzing late-stage reaction: {rsmi}")
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                # Check if we have at least 2 reactants
                if len(reactants) >= 2:
                    print(f"Found {len(reactants)} reactants in late-stage reaction")

                    # Check complexity of reactants
                    complex_fragments = 0
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            # Consider a fragment complex if it has at least 10 atoms
                            # and contains at least one ring
                            if (
                                mol
                                and mol.GetNumAtoms() >= 10
                                and mol.GetRingInfo().NumRings() > 0
                            ):
                                complex_fragments += 1
                                print(
                                    f"Found complex fragment: {reactant} with {mol.GetNumAtoms()} atoms and {mol.GetRingInfo().NumRings()} rings"
                                )
                        except Exception as e:
                            print(f"Error processing reactant {reactant}: {e}")

                    # Check if this is a coupling reaction
                    is_coupling = False
                    coupling_reactions = [
                        "Suzuki",
                        "Negishi",
                        "Stille",
                        "Heck",
                        "Sonogashira",
                        "Buchwald-Hartwig",
                        "Ullmann",
                        "Kumada",
                    ]

                    for rxn_type in coupling_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Detected {rxn_type} coupling reaction")
                            is_coupling = True
                            break

                    # If not a named coupling reaction, check for C-C bond formation
                    if not is_coupling:
                        # Check if reaction involves C-C bond formation between two complex fragments
                        product = rsmi.split(">")[-1]
                        product_mol = Chem.MolFromSmiles(product)

                        if product_mol and complex_fragments >= 2:
                            print(f"Found reaction with at least 2 complex fragments")
                            is_coupling = True

                    if complex_fragments >= 2 and is_coupling:
                        print(
                            "Found convergent synthesis with late-stage coupling of complex fragments"
                        )
                        result = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return result

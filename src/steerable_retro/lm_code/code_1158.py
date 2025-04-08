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
    Detects if the synthesis involves late-stage coupling of major fragments.
    """
    fragment_paths = []
    max_depth = 0

    def dfs_traverse(node, path=None, depth=0):
        nonlocal fragment_paths, max_depth

        if path is None:
            path = []

        if depth > max_depth:
            max_depth = depth

        # Create a copy of the current path
        current_path = path.copy()

        # Add current node to path
        current_path.append((node, depth))

        # If this is a leaf node (starting material)
        if node.get("type") == "mol" and node.get("in_stock", False):
            fragment_paths.append(current_path)
            return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, current_path, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Find coupling reactions
    coupling_reactions = []

    for i in range(len(fragment_paths)):
        for j in range(i + 1, len(fragment_paths)):
            path1 = fragment_paths[i]
            path2 = fragment_paths[j]

            # Find common reaction nodes
            common_reactions = []
            for node1, depth1 in path1:
                for node2, depth2 in path2:
                    if node1 is node2 and node1.get("type") == "reaction":
                        common_reactions.append((node1, depth1))

            if common_reactions:
                # Sort by depth (descending) to find the deepest common reaction
                common_reactions.sort(key=lambda x: x[1], reverse=True)
                deepest_common = common_reactions[0]

                # Check if it's a coupling reaction
                reaction_node = deepest_common[0]
                if "metadata" in reaction_node and "rsmi" in reaction_node["metadata"]:
                    rsmi = reaction_node["metadata"]["rsmi"]

                    # Check if it's a known coupling reaction
                    is_coupling = False
                    coupling_rxn_types = [
                        "Suzuki",
                        "Negishi",
                        "Stille",
                        "Heck",
                        "Sonogashira",
                        "Buchwald-Hartwig",
                        "Ullmann",
                        "N-arylation",
                        "Kumada",
                        "Hiyama-Denmark",
                        "decarboxylative_coupling",
                        "Catellani",
                        "Aryllithium cross-coupling",
                    ]

                    for rxn_type in coupling_rxn_types:
                        if checker.check_reaction(rxn_type, rsmi):
                            is_coupling = True
                            print(f"Found coupling reaction: {rxn_type}")
                            break

                    # If it's not a named coupling, check for C-C bond formation
                    if not is_coupling:
                        # Extract reactants and product
                        try:
                            reactants = rsmi.split(">")[0].split(".")
                            product = rsmi.split(">")[-1]

                            # Check if both reactants are substantial fragments
                            # (at least 8 atoms each is a reasonable threshold)
                            substantial_fragments = True
                            for reactant in reactants:
                                mol = Chem.MolFromSmiles(reactant)
                                if mol and mol.GetNumAtoms() < 8:
                                    substantial_fragments = False

                            if substantial_fragments:
                                is_coupling = True
                                print("Found coupling of substantial fragments")
                        except Exception as e:
                            print(f"Error analyzing reactants: {e}")

                    if is_coupling:
                        coupling_reactions.append(deepest_common)

    # Sort coupling reactions by depth (ascending)
    coupling_reactions.sort(key=lambda x: x[1])

    # Check if we have coupling reactions
    if not coupling_reactions:
        print("No fragment coupling reactions found")
        return False

    # Check if the first coupling reaction is in the second half of the synthesis
    # (lower depth values correspond to later stages in synthesis)
    first_coupling = coupling_reactions[0]
    is_late_stage = first_coupling[1] <= (max_depth / 2)

    print(f"First fragment coupling occurs at depth {first_coupling[1]} (max depth: {max_depth})")
    print(f"Late-stage coupling: {is_late_stage}")

    return is_late_stage

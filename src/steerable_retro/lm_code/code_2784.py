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
    Detects a convergent synthesis strategy where multiple complex fragments
    are combined in late-stage reactions.
    """
    convergent_detected = False

    def dfs_traverse(node, current_depth=0):
        nonlocal convergent_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Extract depth from ID if available, otherwise use current_depth
            node_depth = current_depth
            depth_match = re.search(r"Depth: (\d+)", node.get("metadata", {}).get("ID", ""))
            if depth_match:
                node_depth = int(depth_match.group(1))

            print(f"Examining reaction at depth {node_depth}, with {len(reactants)} reactants")

            # Check if this is a late-stage reaction (depth 0, 1, 2, or 3)
            if node_depth <= 3:
                # Check if we have multiple reactants
                if len(reactants) >= 2:
                    complex_fragments = 0
                    non_complex_fragments = []

                    for r in reactants:
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            # Define "complex" as having significant structure
                            ring_info = mol.GetRingInfo()
                            num_rings = ring_info.NumRings()
                            num_atoms = mol.GetNumHeavyAtoms()

                            # A complex fragment has either multiple rings or is relatively large
                            if (num_rings >= 2 and num_atoms >= 8) or num_atoms >= 12:
                                complex_fragments += 1
                                print(
                                    f"Found complex fragment: {r} with {num_rings} rings and {num_atoms} atoms"
                                )
                            else:
                                non_complex_fragments.append((r, num_rings, num_atoms))

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
                        "Chan-Lam",
                        "Kumada",
                        "Hiyama-Denmark",
                        "N-arylation",
                    ]

                    for rxn_name in coupling_reactions:
                        if checker.check_reaction(rxn_name, rsmi):
                            is_coupling = True
                            print(f"Detected {rxn_name} coupling reaction")
                            break

                    # Consider convergent if:
                    # 1. We have multiple complex fragments regardless of reaction type
                    # 2. We have a coupling reaction with at least one complex fragment
                    # 3. We have one complex fragment and at least one meaningful fragment (≥5 atoms)
                    has_meaningful_fragment = any(
                        atoms >= 5 for _, _, atoms in non_complex_fragments
                    )

                    if complex_fragments >= 2:
                        print(
                            f"Detected convergent synthesis with {complex_fragments} complex fragments"
                        )
                        convergent_detected = True
                    elif is_coupling and complex_fragments >= 1:
                        print(
                            f"Detected convergent synthesis with coupling reaction and {complex_fragments} complex fragment"
                        )
                        convergent_detected = True
                    elif complex_fragments >= 1 and has_meaningful_fragment:
                        print(
                            f"Detected convergent synthesis with {complex_fragments} complex fragment and meaningful smaller fragments"
                        )
                        convergent_detected = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)
    print(f"Final convergent_detected status: {convergent_detected}")
    return convergent_detected

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
    This function detects if the synthesis involves an early multicomponent reaction
    with 3 or more components converging to form a complex intermediate.
    """
    found_multicomponent = False

    def dfs_traverse(node, depth=0):
        nonlocal found_multicomponent

        if found_multicomponent:
            return  # Early return if we already found what we're looking for

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            # Check if this is an early-stage reaction (depth 2+)
            if depth >= 2:
                print(
                    f"Examining reaction at depth {depth} with {len(reactants)} reactants: {rsmi}"
                )

                # Check for known multicomponent reactions first
                if (
                    checker.check_reaction("Ugi reaction", rsmi)
                    or checker.check_reaction("A3 coupling", rsmi)
                    or checker.check_reaction("A3 coupling to imidazoles", rsmi)
                    or checker.check_reaction(
                        "Petasis reaction with amines and boronic acids", rsmi
                    )
                    or checker.check_reaction(
                        "Petasis reaction with amines and boronic esters", rsmi
                    )
                    or checker.check_reaction(
                        "Petasis reaction with amines aldehydes and boronic acids", rsmi
                    )
                ):
                    print(f"Found known multicomponent reaction at depth {depth}")
                    found_multicomponent = True
                    return

                # Count number of distinct reactants (at least 3 for multicomponent)
                if len(reactants) >= 3:
                    print(
                        f"Found potential multicomponent reaction with {len(reactants)} components at depth {depth}"
                    )

                    # Validate reactants are chemically distinct by converting to canonical SMILES
                    unique_reactants = set()
                    for r in reactants:
                        # Skip empty reactants
                        if not r.strip():
                            continue

                        # Remove atom mapping for comparison
                        clean_r = r
                        for match in re.finditer(r":[0-9]+", r):
                            clean_r = clean_r.replace(match.group(0), "")

                        try:
                            mol = Chem.MolFromSmiles(clean_r)
                            if mol:
                                canonical_smiles = Chem.MolToSmiles(mol)
                                unique_reactants.add(canonical_smiles)
                                print(f"  Unique reactant: {canonical_smiles}")
                        except Exception as e:
                            # Skip invalid SMILES
                            print(f"  Error processing reactant {r}: {e}")
                            continue

                    print(f"Found {len(unique_reactants)} chemically distinct reactants")
                    if len(unique_reactants) >= 3:
                        print(
                            f"Confirmed multicomponent reaction with {len(unique_reactants)} chemically distinct reactants"
                        )
                        found_multicomponent = True
                        return

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {found_multicomponent}")
    return found_multicomponent

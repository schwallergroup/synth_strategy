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
    This function detects if the synthesis uses a convergent approach where
    two complex fragments are joined in the final step.
    """
    convergent_final_step = False

    def dfs_traverse(node):
        nonlocal convergent_final_step

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            depth = node.get("metadata", {}).get("depth", None)
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]

            # Check if this is a late-stage step (depth 0, 1, or 2)
            if depth is None or depth in [0, "0", 1, "1", 2, "2"]:
                reactants = reactants_part.split(".")

                # Count complex reactants (those with significant structure)
                complex_reactants = 0
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and (
                            mol.GetNumAtoms() > 10
                            or rdMolDescriptors.CalcNumRings(mol) >= 2
                        ):
                            complex_reactants += 1
                            print(
                                f"Complex reactant found: {reactant} with {mol.GetNumAtoms()} atoms and {rdMolDescriptors.CalcNumRings(mol)} rings"
                            )
                    except Exception as e:
                        print(f"Error processing reactant {reactant}: {e}")
                        continue

                # Check if this is a coupling reaction commonly used in convergent synthesis
                is_coupling_reaction = False
                if rsmi:
                    try:
                        is_coupling_reaction = (
                            checker.check_reaction(
                                "Suzuki coupling with boronic acids", rsmi
                            )
                            or checker.check_reaction(
                                "Suzuki coupling with boronic esters", rsmi
                            )
                            or checker.check_reaction("Negishi coupling", rsmi)
                            or checker.check_reaction("Stille reaction_aryl", rsmi)
                            or checker.check_reaction("Heck terminal vinyl", rsmi)
                            or checker.check_reaction(
                                "Sonogashira alkyne_aryl halide", rsmi
                            )
                        )
                        if is_coupling_reaction:
                            print(
                                f"Detected coupling reaction at depth {depth}: {rsmi}"
                            )
                    except Exception as e:
                        print(f"Error checking reaction type: {e}")

                # Determine if this is a convergent step
                if complex_reactants >= 2:
                    print(
                        f"Detected convergent synthesis with {complex_reactants} complex fragments at depth {depth}"
                    )
                    convergent_final_step = True
                elif complex_reactants >= 1 and is_coupling_reaction:
                    print(
                        f"Detected convergent synthesis with coupling reaction at depth {depth}"
                    )
                    convergent_final_step = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    if not convergent_final_step:
        print("No convergent synthesis pattern detected")

    return convergent_final_step

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
    This function detects a late-stage fragment coupling strategy where two complex
    fragments are joined in the final synthetic step.
    """
    print("Starting late_stage_fragment_coupling_strategy analysis...")
    final_coupling_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal final_coupling_detected

        if node["type"] == "reaction":
            # Extract reaction SMILES
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check if this is a late-stage reaction (depth 0 or 1)
                if depth <= 1:
                    print(f"Analyzing late-stage reaction at depth {depth}: {rsmi}")

                    # Check for coupling reactions
                    coupling_reactions = [
                        "Suzuki",
                        "Negishi",
                        "Stille",
                        "Heck",
                        "Sonogashira",
                        "Buchwald-Hartwig",
                        "Ullmann-Goldberg",
                        "N-arylation",
                        "Hiyama-Denmark Coupling",
                        "Kumada cross-coupling",
                        "Aryllithium cross-coupling",
                        "decarboxylative_coupling",
                    ]

                    for rxn_type in coupling_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Found late-stage coupling reaction: {rxn_type}")
                            final_coupling_detected = True
                            return

                    # If no specific coupling reaction found, check for complex fragments
                    if not final_coupling_detected and len(reactants) >= 2:
                        complex_reactants = 0
                        for reactant in reactants:
                            try:
                                mol = Chem.MolFromSmiles(reactant)
                                if mol:
                                    # Define complex fragment: >12 atoms or ≥2 rings or contains specific functional groups
                                    num_atoms = mol.GetNumAtoms()
                                    num_rings = mol.GetRingInfo().NumRings()

                                    has_complex_fg = (
                                        checker.check_fg("Aromatic halide", reactant)
                                        or checker.check_fg("Boronic acid", reactant)
                                        or checker.check_fg("Boronic ester", reactant)
                                        or checker.check_ring("benzene", reactant)
                                        or checker.check_ring("pyridine", reactant)
                                        or checker.check_ring("indole", reactant)
                                    )

                                    if num_atoms > 12 or num_rings >= 2 or has_complex_fg:
                                        complex_reactants += 1
                                        print(
                                            f"Found complex fragment: {reactant} (atoms: {num_atoms}, rings: {num_rings})"
                                        )
                            except Exception as e:
                                print(f"Error analyzing reactant {reactant}: {e}")

                        if complex_reactants >= 2:
                            print(
                                f"Found late-stage fragment coupling with {complex_reactants} complex fragments"
                            )
                            final_coupling_detected = True
                            return

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Late-stage fragment coupling detected: {final_coupling_detected}")
    return final_coupling_detected

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
    Detects if the synthesis route uses a convergent strategy with
    a final coupling step combining two complex fragments.
    """
    convergent_coupling = False

    def dfs_traverse(node):
        nonlocal convergent_coupling

        if node["type"] == "reaction" and node.get("depth", 0) == 0:  # Check only final step
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check if this is a coupling reaction
                is_coupling = False
                coupling_reactions = [
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic esters",
                    "Buchwald-Hartwig",
                    "Negishi coupling",
                    "Stille reaction_aryl",
                    "Sonogashira alkyne_aryl halide",
                    "Heck terminal vinyl",
                    "Ullmann condensation",
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                ]

                for rxn_type in coupling_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_coupling = True
                        print(f"Detected {rxn_type} as final coupling step")
                        break

                # If not a known coupling reaction, check for C-C, C-N, or C-O bond formation
                if not is_coupling:
                    # Check if there's a new bond formed between fragments
                    # This is a simplified check - in a real implementation,
                    # we would need to analyze the atom mapping more carefully
                    if len(reactants) >= 2:
                        is_coupling = True
                        print("Detected potential coupling reaction in final step")

                if is_coupling:
                    # Check if there are at least 2 complex reactants
                    complex_reactants = 0
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            # Define complexity as having at least one ring
                            ring_info = reactant_mol.GetRingInfo()
                            if ring_info.NumRings() > 0:
                                complex_reactants += 1

                            # Alternative complexity measure: number of atoms
                            elif reactant_mol.GetNumHeavyAtoms() >= 8:
                                complex_reactants += 1

                    if complex_reactants >= 2:
                        convergent_coupling = True
                        print(
                            f"Confirmed convergent coupling with {complex_reactants} complex fragments"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return convergent_coupling

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
    This function detects if the synthetic route employs a convergent synthesis approach
    where multiple complex fragments are combined in late-stage reactions.
    """
    # Track if we found a convergent synthesis pattern
    convergent_synthesis_found = False

    def has_functional_groups(mol_smiles):
        """Check if molecule has significant functional groups"""
        functional_groups = [
            "Carboxylic acid",
            "Ester",
            "Amide",
            "Amine",
            "Alcohol",
            "Aldehyde",
            "Ketone",
            "Nitrile",
            "Halide",
            "Boronic acid",
        ]
        for fg in functional_groups:
            if checker.check_fg(fg, mol_smiles):
                return True
        return False

    def is_coupling_reaction(rxn_smiles):
        """Check if reaction is a coupling reaction typically used in convergent synthesis"""
        coupling_reactions = [
            "Suzuki",
            "Negishi",
            "Stille",
            "Heck",
            "Sonogashira",
            "Buchwald-Hartwig",
            "Ullmann",
            "Wittig",
            "Grignard",
            "Diels-Alder",
        ]
        for rxn in coupling_reactions:
            if checker.check_reaction(rxn, rxn_smiles):
                print(f"Detected {rxn} coupling reaction")
                return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_synthesis_found

        # Store depth in metadata for reference
        if "metadata" not in node:
            node["metadata"] = {}
        node["metadata"]["depth"] = depth

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if we're combining multiple complex fragments
            complex_fragments = 0
            complex_reactants = []

            for reactant in reactants:
                if not reactant:
                    continue

                mol = Chem.MolFromSmiles(reactant)
                if not mol:
                    continue

                # Define complexity based on size, rings, and functional groups
                is_complex = mol.GetNumAtoms() > 12 and (
                    mol.GetRingInfo().NumRings() > 0 or has_functional_groups(reactant)
                )

                if is_complex:
                    complex_fragments += 1
                    complex_reactants.append(reactant)
                    print(
                        f"Complex fragment found: {reactant} (atoms: {mol.GetNumAtoms()}, rings: {mol.GetRingInfo().NumRings()})"
                    )

            # Criteria for convergent synthesis:
            # 1. Combining 2+ complex fragments
            # 2. In a late-stage reaction (depth ≤ 3)
            # 3. Using a coupling reaction (optional but strengthens evidence)
            is_late_stage = depth <= 3
            is_coupling = is_coupling_reaction(rsmi)

            if complex_fragments >= 2 and is_late_stage:
                print(
                    f"Convergent synthesis detected: combining {complex_fragments} complex fragments in {'late' if is_late_stage else 'early'}-stage reaction"
                )
                print(f"Complex reactants: {complex_reactants}")
                if is_coupling:
                    print("Using a coupling reaction, strong evidence for convergent synthesis")
                convergent_synthesis_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Convergent synthesis found: {convergent_synthesis_found}")
    return convergent_synthesis_found

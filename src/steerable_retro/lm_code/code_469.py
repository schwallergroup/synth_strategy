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
    Detects if the synthesis follows a linear strategy where each step builds upon
    a single precursor molecule rather than combining multiple complex fragments.
    """
    is_linear = True

    def is_simple_reagent(mol_smiles):
        """Identify common simple reagents that shouldn't count as complex reactants"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return True

        # Small molecules are considered simple reagents
        if mol.GetNumHeavyAtoms() <= 6:
            return True

        # Check for common reagent functional groups
        common_reagent_fgs = [
            "Triflate",
            "Tosylate",
            "Mesylate",
            "Primary halide",
            "Secondary halide",
            "Tertiary halide",
            "Acyl halide",
            "Boronic acid",
            "Boronic ester",
            "Zinc halide",
            "Tin",
            "Alkyl lithium",
            "Aryl lithium",
            "Magnesium halide",
        ]

        for fg in common_reagent_fgs:
            if checker.check_fg(fg, mol_smiles):
                return True

        # Check for common simple reagent rings
        common_reagent_rings = ["cyclopropane", "cyclobutane", "oxirane", "aziridine", "thiirane"]
        for ring in common_reagent_rings:
            if checker.check_ring(ring, mol_smiles) and mol.GetNumHeavyAtoms() <= 10:
                return True

        return False

    def is_coupling_reaction(rsmi):
        """Check if this is a coupling reaction that inherently needs two components"""
        coupling_rxn_types = [
            "Suzuki coupling",
            "Buchwald-Hartwig",
            "Sonogashira",
            "Heck",
            "Stille reaction",
            "Negishi coupling",
            "Ullmann condensation",
            "Kumada cross-coupling",
            "Hiyama-Denmark Coupling",
            "Chan-Lam",
            "Catellani reaction",
            "Goldberg coupling",
            "Ullmann-Goldberg Substitution",
        ]

        for rxn_type in coupling_rxn_types:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"Detected coupling reaction: {rxn_type}")
                return True

        return False

    def is_acceptable_multicomponent_reaction(rsmi):
        """Check if this is a multi-component reaction that is commonly used in linear synthesis"""
        acceptable_rxn_types = [
            "Reductive amination",
            "Wittig",
            "Aldol condensation",
            "Schotten-Baumann",
            "Ugi reaction",
            "Petasis reaction",
            "Pictet-Spengler",
            "Michael addition",
            "Knoevenagel Condensation",
            "Henry Reaction",
            "Mitsunobu",
            "Eschweiler-Clarke",
            "Paal-Knorr pyrrole synthesis",
        ]

        for rxn_type in acceptable_rxn_types:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"Detected acceptable multi-component reaction: {rxn_type}")
                return True

        return False

    def dfs_traverse(node):
        nonlocal is_linear

        if not is_linear:
            return  # Stop traversal if already determined to be non-linear

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Count non-trivial reactants (excluding simple reagents)
            complex_reactants = []
            for reactant_smiles in reactants_smiles:
                try:
                    if not is_simple_reagent(reactant_smiles):
                        mol = Chem.MolFromSmiles(reactant_smiles)
                        if mol and mol.GetNumHeavyAtoms() > 8:  # Threshold for complex molecules
                            complex_reactants.append(reactant_smiles)
                except Exception as e:
                    print(f"Error processing reactant {reactant_smiles}: {e}")

            # If more than one complex reactant, check if it's an acceptable reaction type
            if len(complex_reactants) > 1:
                if not (is_coupling_reaction(rsmi) or is_acceptable_multicomponent_reaction(rsmi)):
                    print(
                        f"Non-linear step detected with {len(complex_reactants)} complex reactants"
                    )
                    is_linear = False
                    return
                else:
                    print("Multi-component reaction is acceptable in linear synthesis")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return is_linear

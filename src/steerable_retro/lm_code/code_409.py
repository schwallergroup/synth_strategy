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
    This function detects a convergent synthesis strategy where two or more
    complex fragments are joined in the final steps of the synthesis.
    """
    is_convergent = False

    def calculate_depth(node):
        """Calculate depth if not provided in metadata"""
        max_depth = 0
        for child in node.get("children", []):
            if child["type"] == "reaction":
                child_depth = child["metadata"].get("depth", calculate_depth(child))
                max_depth = max(max_depth, child_depth + 1)
            elif child["type"] == "mol" and not child.get("in_stock", False):
                for grandchild in child.get("children", []):
                    if grandchild["type"] == "reaction":
                        grandchild_depth = grandchild["metadata"].get(
                            "depth", calculate_depth(grandchild)
                        )
                        max_depth = max(max_depth, grandchild_depth + 1)
        return max_depth

    def has_complex_structure(mol):
        """Check if molecule has complex structural features"""
        mol_smiles = Chem.MolToSmiles(mol)

        # Check for rings
        ring_count = 0
        for ring_name in [
            "benzene",
            "pyridine",
            "piperidine",
            "cyclohexane",
            "furan",
            "pyrrole",
            "thiophene",
            "morpholine",
            "indole",
            "pyrazole",
        ]:
            if checker.check_ring(ring_name, mol_smiles):
                ring_count += 1
                print(f"Found ring: {ring_name}")
                if ring_count >= 1:  # Even one ring can indicate complexity
                    return True

        # Check for functional groups that indicate complexity
        complex_fg_count = 0
        for fg_name in [
            "Ester",
            "Amide",
            "Ether",
            "Alcohol",
            "Amine",
            "Nitrile",
            "Carboxylic acid",
            "Aromatic halide",
            "Primary halide",
            "Secondary halide",
            "Tertiary halide",
        ]:
            if checker.check_fg(fg_name, mol_smiles):
                complex_fg_count += 1
                print(f"Found functional group: {fg_name}")
                if complex_fg_count >= 1:  # Even one complex FG can indicate complexity
                    return True

        return False

    def is_coupling_reaction(rsmi):
        """Check if the reaction is a coupling reaction"""
        coupling_reactions = [
            "Suzuki coupling with boronic acids",
            "Suzuki coupling with boronic esters",
            "Negishi coupling",
            "Stille reaction",
            "Heck terminal vinyl",
            "Sonogashira acetylene_aryl halide",
            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
            "Ullmann condensation",
            "Williamson Ether Synthesis",  # Added this important coupling reaction
            "S-alkylation of thiols",
            "N-alkylation of primary amines with alkyl halides",
            "N-alkylation of secondary amines with alkyl halides",
        ]

        for rxn_name in coupling_reactions:
            if checker.check_reaction(rxn_name, rsmi):
                print(f"Detected coupling reaction: {rxn_name}")
                return True

        # Check for general etherification (which includes Williamson)
        if "O" in rsmi and ("Br" in rsmi or "Cl" in rsmi or "I" in rsmi):
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product has an ether group
            if checker.check_fg("Ether", product):
                print("Detected potential etherification reaction")
                return True

        return False

    def dfs_traverse(node, current_depth=0):
        nonlocal is_convergent

        if node["type"] == "reaction":
            # Get or calculate depth
            depth = node["metadata"].get("depth", current_depth)

            # Focus on late-stage reactions (depth 0, 1, or 2)
            if depth <= 2:
                # Extract reaction SMILES
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    return

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Split into reactants and product
                parts = rsmi.split(">")
                if len(parts) < 3:
                    return

                reactants = parts[0].split(".")

                # Count complex reactants (those with significant structure)
                complex_reactants = 0
                total_atom_count = 0

                for r in reactants:
                    if r:
                        try:
                            mol = Chem.MolFromSmiles(r)
                            if mol:
                                atom_count = mol.GetNumAtoms()
                                total_atom_count += atom_count

                                # Consider a molecule complex if it has >10 atoms or complex structure
                                is_complex = atom_count > 10 or has_complex_structure(mol)
                                print(f"Reactant: {r}, Atoms: {atom_count}, Complex: {is_complex}")
                                if is_complex:
                                    complex_reactants += 1
                        except Exception as e:
                            print(f"Error processing reactant {r}: {e}")

                # Check if this is a coupling reaction
                is_coupling = is_coupling_reaction(rsmi)

                # Criteria for convergent synthesis:
                # 1. At least 2 complex reactants
                # 2. A coupling reaction with at least 1 complex reactant
                # 3. A late-stage reaction (depth ≤ 1) with high total atom count
                if (
                    complex_reactants >= 2
                    or (is_coupling and complex_reactants >= 1)
                    or (depth <= 1 and total_atom_count > 30)
                ):
                    print(
                        f"Detected convergent synthesis at depth {depth} with {complex_reactants} complex fragments"
                    )
                    print(f"Total atom count: {total_atom_count}")
                    is_convergent = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: {is_convergent}")
    return is_convergent

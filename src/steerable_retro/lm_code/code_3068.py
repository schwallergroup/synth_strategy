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
    Detects a synthesis strategy that is linear until the final step,
    where it becomes convergent by combining two complex fragments.
    """
    # Initialize tracking variables
    linear_steps = []
    convergent_steps = []
    main_pathway = []

    def is_simple_reagent(smiles):
        """Check if a molecule is a simple reagent like solvent, catalyst, or common reagent"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        # Check for common functional groups that indicate reagents
        common_reagent_fgs = [
            "Acyl halide",
            "Triflate",
            "Mesylate",
            "Tosylate",
            "Magnesium halide",
            "Zinc halide",
            "Alkyl lithium",
            "Tin",
            "Silane",
            "Boronic acid",
            "Boronic ester",
            "Grignard",
            "Diazo",
        ]

        for fg in common_reagent_fgs:
            if checker.check_fg(fg, smiles):
                return True

        # Simple reagents typically have few atoms
        if mol.GetNumAtoms() <= 7:
            return True

        # Check for common simple reagents by pattern
        common_patterns = [
            # Solvents
            "O",
            "CCO",
            "CCCO",
            "CCCCO",
            "CC(C)O",
            "C1CCCCC1",
            "CC(=O)C",
            "CN(C)C=O",
            "CS(=O)C",
            "ClCCl",
            "ClCCCl",
            "CC#N",
            "C=O",
            # Gases
            "[H][H]",
            "C(=O)=O",
            "N#N",
            "O=O",
            "C",
            # Acids/Bases
            "CC(=O)O",
            "O=C(O)C(F)(F)F",
            "C1=CC=C(C=C1)S(=O)(=O)O",
            "CN1C=NC=C1",
            "CC1=CN=C[NH]1",
            "CN1C=CN=C1",
            "C1=CC=NC=C1",
            # Common reagents
            "Br",
            "Cl",
            "I",
            "F",
            "[Na+]",
            "[K+]",
            "[Li+]",
            "B(O)O",
            "P(Cl)Cl",
        ]

        for pattern in common_patterns:
            if Chem.MolFromSmiles(pattern) and Chem.MolFromSmiles(pattern).GetSubstructMatch(mol):
                return True

        return False

    def is_complex_fragment(smiles, node=None):
        """Check if a molecule is a complex fragment (not a simple reagent and structurally significant)"""
        if is_simple_reagent(smiles):
            return False

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        # Complex fragments typically have more atoms
        if mol.GetNumAtoms() < 8:
            return False

        # Check if it's in-stock and not complex enough
        if node and node.get("in_stock", False) and mol.GetNumAtoms() < 12:
            return False

        # Check for presence of rings which often indicate complexity
        ring_count = 0
        for ring in [
            "benzene",
            "pyridine",
            "piperidine",
            "cyclohexane",
            "cyclopentane",
            "furan",
            "pyrrole",
        ]:
            if checker.check_ring(ring, smiles):
                ring_count += 1

        if ring_count > 0:
            return True

        # Check for multiple functional groups which indicate complexity
        fg_count = 0
        for fg in [
            "Ester",
            "Amide",
            "Carboxylic acid",
            "Alcohol",
            "Amine",
            "Nitrile",
            "Nitro group",
        ]:
            if checker.check_fg(fg, smiles):
                fg_count += 1

        if fg_count >= 2:
            return True

        # Default to atom count for other cases
        return mol.GetNumAtoms() > 15

    def get_mol_node_from_smiles(node, smiles):
        """Find a molecule node in children by its SMILES"""
        for child in node.get("children", []):
            if child["type"] == "mol" and child["smiles"] == smiles:
                return child
        return None

    def dfs_traverse(node, depth=0, path=None):
        if path is None:
            path = []

        current_path = path + [node]

        # If this is a reaction node, analyze it
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Count complex reactants
                complex_reactants = []
                for r_smiles in reactants_smiles:
                    r_node = get_mol_node_from_smiles(node, r_smiles)
                    if is_complex_fragment(r_smiles, r_node):
                        complex_reactants.append(r_smiles)

                # Determine if this is a linear or convergent step
                if len(complex_reactants) == 1:
                    linear_steps.append((depth, node, complex_reactants[0]))
                    print(
                        f"Found linear step at depth {depth} with complex reactant: {complex_reactants[0][:30]}..."
                    )
                elif len(complex_reactants) >= 2:
                    convergent_steps.append((depth, node, complex_reactants))
                    print(
                        f"Found convergent step at depth {depth} with {len(complex_reactants)} complex reactants"
                    )

                # Track the main synthetic pathway (the one with the most complex steps)
                if len(complex_reactants) > 0:
                    main_pathway.append((depth, node))
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

    # Start traversal
    dfs_traverse(route)

    # Sort steps by depth
    linear_steps.sort(key=lambda x: x[0])
    convergent_steps.sort(key=lambda x: x[0])
    main_pathway.sort(key=lambda x: x[0])

    # Check if we have a linear-to-convergent pattern
    has_linear_to_convergent = False

    if linear_steps and convergent_steps:
        # The final step should be convergent (lowest depth)
        final_step = min(main_pathway, key=lambda x: x[0]) if main_pathway else None

        if final_step and final_step in [step[:2] for step in convergent_steps]:
            # Check if we have at least 2 linear steps before the convergent step
            linear_before_convergent = [step for step in linear_steps if step[0] > final_step[0]]

            if len(linear_before_convergent) >= 2:
                has_linear_to_convergent = True
                print(
                    f"Found linear-to-convergent pattern with {len(linear_before_convergent)} linear steps before final convergent step"
                )

    print(f"Linear steps: {len(linear_steps)}, Convergent steps: {len(convergent_steps)}")
    print(f"Linear-to-convergent pattern: {has_linear_to_convergent}")

    return has_linear_to_convergent

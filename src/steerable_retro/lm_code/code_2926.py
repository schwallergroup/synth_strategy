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
    This function detects if the synthesis follows a convergent approach where
    two complex fragments (>10 atoms) are combined in a late-stage reaction.
    """
    convergent_found = False

    def is_coupling_reaction(rsmi):
        """Check if the reaction is a coupling reaction type"""
        # Common coupling reactions
        coupling_reactions = [
            "Suzuki coupling with boronic acids",
            "Suzuki coupling with boronic esters",
            "Amidation",
            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
            "Esterification of Carboxylic Acids",
            "Williamson Ether Synthesis",
            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
            "Heck terminal vinyl",
            "Negishi coupling",
            "Stille reaction_aryl",
            "Sonogashira alkyne_aryl halide",
        ]

        for reaction_type in coupling_reactions:
            if checker.check_reaction(reaction_type, rsmi):
                print(f"Detected coupling reaction: {reaction_type}")
                return True

        # If no specific coupling reaction is found, check for general bond formation
        try:
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants_mol = Chem.MolFromSmiles(reactants_part)
            product_mol = Chem.MolFromSmiles(product_part)

            if reactants_mol and product_mol:
                # If product has fewer molecules than reactants, it's likely a coupling
                reactant_count = len(reactants_part.split("."))
                if "." not in product_part and reactant_count >= 2:
                    print("Detected general coupling: multiple reactants forming single product")
                    return True
        except Exception as e:
            print(f"Error checking for general coupling: {e}")

        return False

    def contributes_significantly(reactant, rsmi):
        """Check if the reactant contributes significantly to the product structure"""
        try:
            product = rsmi.split(">")[-1]
            reactant_mol = Chem.MolFromSmiles(reactant)
            product_mol = Chem.MolFromSmiles(product)

            if reactant_mol and product_mol:
                reactant_atoms = reactant_mol.GetNumHeavyAtoms()
                product_atoms = product_mol.GetNumHeavyAtoms()

                # Reactant should contribute at least 30% of product's atoms
                contribution = reactant_atoms / product_atoms
                print(
                    f"Reactant contribution: {contribution:.2f} ({reactant_atoms}/{product_atoms})"
                )
                return contribution >= 0.3
        except Exception as e:
            print(f"Error checking reactant contribution: {e}")

        return False

    def dfs_traverse(node, current_depth=0):
        nonlocal convergent_found

        if node["type"] == "reaction" and "metadata" in node:
            # Extract depth safely
            depth = current_depth
            if "ID" in node["metadata"] and "Depth:" in node["metadata"]["ID"]:
                try:
                    depth_str = node["metadata"]["ID"].split("Depth: ")[-1]
                    # Handle case where there might be text after the depth number
                    depth = int("".join(c for c in depth_str if c.isdigit()))
                except (ValueError, IndexError):
                    print(
                        f"Could not parse depth from ID: {node['metadata'].get('ID')}, using calculated depth: {current_depth}"
                    )

            # Check if it's a late-stage reaction (depth ≤ 3)
            if depth <= 3:
                print(f"Examining late-stage reaction at depth {depth}")
                if "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    reactants_part = rsmi.split(">")[0]
                    reactants = reactants_part.split(".")

                    # Check if this is a coupling reaction
                    if is_coupling_reaction(rsmi):
                        # Count complex fragments (>10 atoms) that contribute significantly
                        complex_fragments = 0
                        for reactant in reactants:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                # Count heavy atoms (non-hydrogen)
                                atom_count = mol.GetNumHeavyAtoms()
                                if atom_count > 10 and contributes_significantly(reactant, rsmi):
                                    complex_fragments += 1
                                    print(
                                        f"Found complex fragment with {atom_count} atoms: {reactant}"
                                    )

                        if complex_fragments >= 2:
                            convergent_found = True
                            print(
                                f"Found convergent synthesis at depth {depth} with {complex_fragments} complex fragments"
                            )
                else:
                    print(f"Reaction node at depth {depth} missing rsmi field")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)
    return convergent_found

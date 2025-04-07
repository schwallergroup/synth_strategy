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
    This function detects if a synthetic route involves multiple functionalization
    steps on a pyrazole scaffold.
    """
    pyrazole_modifications = 0
    pyrazole_present = False

    # Track which molecules contain pyrazole to avoid double-counting
    processed_reactions = set()

    def dfs_traverse(node, depth=0):
        nonlocal pyrazole_modifications, pyrazole_present

        if node["type"] == "mol":
            # Check if this molecule contains a pyrazole ring
            if checker.check_ring("pyrazole", node["smiles"]):
                pyrazole_present = True
                print(f"Depth {depth}: Pyrazole detected in molecule: {node['smiles']}")

        elif node["type"] == "reaction":
            # Skip if we've already processed this reaction
            reaction_id = node["metadata"].get("reaction_hash", "")
            if reaction_id in processed_reactions:
                return
            processed_reactions.add(reaction_id)

            # Extract reactants and product from reaction SMILES
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                # In retrosynthesis, product is the starting material and reactants are the precursors
                if checker.check_ring("pyrazole", product):
                    pyrazole_present = True
                    print(f"Depth {depth}: Pyrazole detected in product: {product}")

                    # Check if any reactant contains pyrazole
                    pyrazole_in_reactants = any(
                        checker.check_ring("pyrazole", r) for r in reactants
                    )

                    # If pyrazole is in both product and at least one reactant, it's being modified
                    if pyrazole_in_reactants:
                        # Check for common functionalization reactions
                        if (
                            checker.check_reaction("Friedel-Crafts acylation", rsmi)
                            or checker.check_reaction("Friedel-Crafts alkylation", rsmi)
                            or checker.check_reaction(
                                "N-alkylation of primary amines with alkyl halides", rsmi
                            )
                            or checker.check_reaction(
                                "N-alkylation of secondary amines with alkyl halides", rsmi
                            )
                            or checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                            or checker.check_reaction(
                                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                            )
                            or checker.check_reaction(
                                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                                rsmi,
                            )
                            or checker.check_reaction("Heck terminal vinyl", rsmi)
                            or checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi)
                            or checker.check_reaction(
                                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                                rsmi,
                            )
                            or checker.check_reaction("Acylation of primary amines", rsmi)
                            or checker.check_reaction("Acylation of secondary amines", rsmi)
                            or checker.check_reaction("Aromatic iodination", rsmi)
                            or checker.check_reaction("Aromatic bromination", rsmi)
                            or checker.check_reaction("Aromatic chlorination", rsmi)
                            or checker.check_reaction("Aromatic fluorination", rsmi)
                        ):

                            pyrazole_modifications += 1
                            print(
                                f"Depth {depth}: Pyrazole modification detected via known reaction. Total modifications: {pyrazole_modifications}"
                            )

                        # Check for functional group changes on pyrazole
                        elif any(
                            checker.check_fg(fg, product)
                            for fg in ["Primary amide", "Secondary amide", "Tertiary amide"]
                        ) and not any(
                            checker.check_fg(fg, r)
                            for fg in ["Primary amide", "Secondary amide", "Tertiary amide"]
                            for r in reactants
                        ):
                            pyrazole_modifications += 1
                            print(
                                f"Depth {depth}: Amide formation on pyrazole detected. Total modifications: {pyrazole_modifications}"
                            )

                        # Check for halogenation
                        elif checker.check_fg("Aromatic halide", product) and not any(
                            checker.check_fg("Aromatic halide", r) for r in reactants
                        ):
                            pyrazole_modifications += 1
                            print(
                                f"Depth {depth}: Halogenation of pyrazole detected. Total modifications: {pyrazole_modifications}"
                            )

                        # Check for iodination specifically
                        elif "I" in product and not any("I" in r for r in reactants):
                            pyrazole_modifications += 1
                            print(
                                f"Depth {depth}: Iodination of pyrazole detected. Total modifications: {pyrazole_modifications}"
                            )

                        # Check for acetylation/acylation specifically
                        elif "C(=O)" in product and not any("C(=O)" in r for r in reactants):
                            pyrazole_modifications += 1
                            print(
                                f"Depth {depth}: Acylation of pyrazole detected. Total modifications: {pyrazole_modifications}"
                            )

                        # Check for Sonogashira coupling (alkyne addition)
                        elif "C#C" in product and not any("C#C" in r for r in reactants):
                            pyrazole_modifications += 1
                            print(
                                f"Depth {depth}: Alkyne addition to pyrazole detected. Total modifications: {pyrazole_modifications}"
                            )

                        # Generic check for any modification based on atom mapping
                        elif ">" in rsmi and "[:" in rsmi:
                            # If we have atom mapping, check if the reaction modifies the pyrazole
                            # This is a fallback for reactions not covered by specific checks
                            pyrazole_modifications += 1
                            print(
                                f"Depth {depth}: Generic pyrazole modification detected. Total modifications: {pyrazole_modifications}"
                            )

                    # Check for pyrazole formation
                    elif not pyrazole_in_reactants:
                        pyrazole_modifications += 1
                        print(
                            f"Depth {depth}: Pyrazole formation detected. Total modifications: {pyrazole_modifications}"
                        )

                # Check for reactions where pyrazole is in reactants but not in product (transformation)
                elif any(
                    checker.check_ring("pyrazole", r) for r in reactants
                ) and not checker.check_ring("pyrazole", product):
                    pyrazole_present = True
                    pyrazole_modifications += 1
                    print(
                        f"Depth {depth}: Pyrazole transformation detected. Total modifications: {pyrazole_modifications}"
                    )

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    print(
        f"Final results: Pyrazole present: {pyrazole_present}, Modifications: {pyrazole_modifications}"
    )
    # Return True if we found pyrazole and multiple modifications to it
    return pyrazole_present and pyrazole_modifications >= 2

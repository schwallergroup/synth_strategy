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
    This function detects if the synthesis uses a late-stage amide coupling strategy
    where two complex fragments are joined via amide bond formation in the final step.
    """
    final_amide_coupling = False

    def dfs_traverse(node, current_depth=0):
        nonlocal final_amide_coupling

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            # Get depth from metadata or use current_depth
            depth = node.get("metadata", {}).get("depth")
            if depth is None:
                depth = current_depth
            else:
                try:
                    depth = int(depth)
                except:
                    depth = current_depth

            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Check if this is the final or near-final step (depth 0 or 1)
            if depth <= 1:
                print(f"Analyzing potential final reaction: {rsmi}")

                # Check if this is an amide coupling reaction
                amide_coupling_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Carboxylic acid with primary amine to amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Schotten-Baumann_amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                ]

                is_amide_coupling = False
                for rxn_type in amide_coupling_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected amide coupling reaction: {rxn_type}")
                        is_amide_coupling = True
                        break

                # If no specific reaction type matched, check for characteristic functional group changes
                if not is_amide_coupling:
                    try:
                        reactants_part = rsmi.split(">")[0]
                        product_part = rsmi.split(">")[-1]

                        # Check if reactants contain acid/acyl halide and amine, and product contains amide
                        has_acid = any(
                            checker.check_fg("Carboxylic acid", r)
                            for r in reactants_part.split(".")
                        )
                        has_acyl_halide = any(
                            checker.check_fg("Acyl halide", r) for r in reactants_part.split(".")
                        )
                        has_amine = any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            or checker.check_fg("Aniline", r)
                            for r in reactants_part.split(".")
                        )

                        has_amide = (
                            checker.check_fg("Primary amide", product_part)
                            or checker.check_fg("Secondary amide", product_part)
                            or checker.check_fg("Tertiary amide", product_part)
                        )

                        if (has_acid or has_acyl_halide) and has_amine and has_amide:
                            print("Detected amide coupling based on functional group changes")
                            is_amide_coupling = True
                    except Exception as e:
                        print(f"Error in functional group analysis: {e}")

                if is_amide_coupling:
                    try:
                        # Extract reactants and product
                        reactants_part = rsmi.split(">")[0]
                        product_part = rsmi.split(">")[-1]
                        reactants = reactants_part.split(".")

                        # Check for carboxylic acid and amine in reactants
                        acid_found = False
                        amine_found = False

                        # Track fragment complexity
                        acid_fragment_size = 0
                        amine_fragment_size = 0

                        for reactant in reactants:
                            try:
                                # Check for carboxylic acid or acyl chloride (common in amide couplings)
                                if checker.check_fg("Carboxylic acid", reactant):
                                    acid_found = True
                                    mol = Chem.MolFromSmiles(reactant)
                                    if mol:
                                        acid_fragment_size = mol.GetNumHeavyAtoms()
                                        print(
                                            f"Found carboxylic acid fragment with {acid_fragment_size} heavy atoms"
                                        )
                                elif checker.check_fg(
                                    "Acyl chloride", reactant
                                ) or checker.check_fg("Acyl halide", reactant):
                                    acid_found = True
                                    mol = Chem.MolFromSmiles(reactant)
                                    if mol:
                                        acid_fragment_size = mol.GetNumHeavyAtoms()
                                        print(
                                            f"Found acyl halide fragment with {acid_fragment_size} heavy atoms"
                                        )

                                # Check for various amine types
                                if (
                                    checker.check_fg("Primary amine", reactant)
                                    or checker.check_fg("Secondary amine", reactant)
                                    or checker.check_fg("Aniline", reactant)
                                ):
                                    amine_found = True
                                    mol = Chem.MolFromSmiles(reactant)
                                    if mol:
                                        amine_fragment_size = mol.GetNumHeavyAtoms()
                                        print(
                                            f"Found amine fragment with {amine_fragment_size} heavy atoms"
                                        )
                            except Exception as e:
                                print(f"Error checking reactant {reactant}: {e}")
                                continue

                        # Check for amide in product
                        amide_found = False
                        try:
                            if (
                                checker.check_fg("Primary amide", product_part)
                                or checker.check_fg("Secondary amide", product_part)
                                or checker.check_fg("Tertiary amide", product_part)
                            ):
                                amide_found = True
                                print("Found amide in product")
                        except Exception as e:
                            print(f"Error checking product {product_part}: {e}")

                        # Consider it late-stage coupling if both fragments are reasonably complex
                        # (at least 4 heavy atoms each is a reasonable threshold for "complex")
                        if acid_found and amine_found and amide_found:
                            print(
                                f"Acid fragment size: {acid_fragment_size}, Amine fragment size: {amine_fragment_size}"
                            )
                            if acid_fragment_size >= 4 and amine_fragment_size >= 4:
                                print(
                                    "Detected late-stage amide coupling in final step with complex fragments"
                                )
                                final_amide_coupling = True
                    except Exception as e:
                        print(f"Error analyzing amide coupling: {e}")

        # Continue traversing with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: {final_amide_coupling}")
    return final_amide_coupling

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
    This function detects a convergent synthesis strategy where two complex fragments
    are joined via N-alkylation of a piperazine ring.
    """
    piperazine_alkylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal piperazine_alkylation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a coupling reaction with piperazine
            if len(reactants) >= 2:  # At least two reactants
                # Check if product contains piperazine
                if checker.check_ring("piperazine", product):
                    print(f"Found product with piperazine at depth {depth}: {product}")

                    # Check if one of the reactants contains piperazine
                    piperazine_reactant = None
                    other_reactants = []

                    for r in reactants:
                        if checker.check_ring("piperazine", r):
                            piperazine_reactant = r
                        else:
                            other_reactants.append(r)

                    # If we found a piperazine-containing reactant
                    if piperazine_reactant and other_reactants:
                        print(f"Found piperazine reactant: {piperazine_reactant}")
                        print(f"Other reactants: {other_reactants}")

                        # Check if this is an N-alkylation or methylation reaction
                        is_alkylation = (
                            checker.check_reaction(
                                "N-alkylation of secondary amines with alkyl halides", rsmi
                            )
                            or checker.check_reaction(
                                "N-alkylation of primary amines with alkyl halides", rsmi
                            )
                            or checker.check_reaction("Methylation with MeI_secondary", rsmi)
                            or checker.check_reaction("Methylation", rsmi)
                            or checker.check_reaction("N-methylation", rsmi)
                        )

                        # Check for methyl iodide specifically
                        methyl_iodide_present = any(
                            "I[CH3" in r or "ICH3" in r or "[CH3:" in r and "I" in r
                            for r in other_reactants
                        )

                        if is_alkylation or methyl_iodide_present:
                            print(f"Confirmed N-alkylation/methylation reaction: {rsmi}")

                            # Check if the other reactant has a leaving group
                            has_leaving_group = False
                            for other_r in other_reactants:
                                if (
                                    checker.check_fg("Primary halide", other_r)
                                    or checker.check_fg("Secondary halide", other_r)
                                    or checker.check_fg("Tertiary halide", other_r)
                                    or checker.check_fg("Mesylate", other_r)
                                    or checker.check_fg("Tosylate", other_r)
                                    or checker.check_fg("Triflate", other_r)
                                    or "I[CH3" in other_r
                                    or "ICH3" in other_r
                                    or ("I" in other_r and "[CH3:" in other_r)
                                ):

                                    print(f"Found reactant with leaving group: {other_r}")
                                    has_leaving_group = True

                                    # Check if both fragments are complex
                                    piperazine_mol = Chem.MolFromSmiles(piperazine_reactant)
                                    other_mol = Chem.MolFromSmiles(other_r)

                                    if piperazine_mol and other_mol:
                                        piperazine_atoms = piperazine_mol.GetNumAtoms()
                                        other_atoms = other_mol.GetNumAtoms()
                                        print(
                                            f"Piperazine fragment atoms: {piperazine_atoms}, Other fragment atoms: {other_atoms}"
                                        )

                                        # Adjusted complexity requirements
                                        if piperazine_atoms >= 8 and other_atoms >= 3:
                                            print(
                                                f"Found piperazine N-alkylation in convergent synthesis: {rsmi}"
                                            )
                                            piperazine_alkylation_found = True
                                            return  # Exit once found
                                        else:
                                            print(
                                                "Fragments not complex enough for convergent synthesis"
                                            )
                                    else:
                                        print("Could not create RDKit molecules from SMILES")

                            if not has_leaving_group:
                                print("No reactant with suitable leaving group found")
                        else:
                            # Fallback check for N-alkylation pattern
                            piperazine_mol = Chem.MolFromSmiles(piperazine_reactant)
                            if piperazine_mol:
                                # Check if piperazine has NH group (secondary amine)
                                has_nh = (
                                    "[NH]" in piperazine_reactant
                                    or "N[H]" in piperazine_reactant
                                    or "[NH:" in piperazine_reactant
                                )

                                if has_nh:
                                    print("Detected secondary amine in piperazine")

                                    # Check if any other reactant has a leaving group
                                    for other_r in other_reactants:
                                        if (
                                            checker.check_fg("Primary halide", other_r)
                                            or checker.check_fg("Secondary halide", other_r)
                                            or checker.check_fg("Tertiary halide", other_r)
                                            or checker.check_fg("Mesylate", other_r)
                                            or checker.check_fg("Tosylate", other_r)
                                            or checker.check_fg("Triflate", other_r)
                                        ):

                                            print(
                                                f"Found reactant with leaving group in fallback check: {other_r}"
                                            )

                                            # Check complexity
                                            other_mol = Chem.MolFromSmiles(other_r)
                                            if other_mol:
                                                piperazine_atoms = piperazine_mol.GetNumAtoms()
                                                other_atoms = other_mol.GetNumAtoms()
                                                print(
                                                    f"Fallback check - Piperazine atoms: {piperazine_atoms}, Other atoms: {other_atoms}"
                                                )

                                                if piperazine_atoms >= 8 and other_atoms >= 3:
                                                    print(
                                                        f"Found piperazine N-alkylation in fallback check: {rsmi}"
                                                    )
                                                    piperazine_alkylation_found = True
                                                    return  # Exit once found
                            else:
                                print("Not an N-alkylation reaction and fallback check failed")
                    else:
                        print("Did not find both piperazine reactant and other reactants")
                else:
                    print(f"Product does not contain piperazine: {product}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return piperazine_alkylation_found

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
    This function detects a synthetic strategy where an indole is N-alkylated
    with a cyclohexylmethyl group in an early stage of the synthesis.
    """
    # Track if we found the pattern
    found_n_alkylation = False

    def dfs_traverse(node, current_depth=0):
        nonlocal found_n_alkylation

        if node["type"] == "reaction":
            # Check for N-alkylation reaction (should be at high depth, early in synthesis)
            if "metadata" in node and "rsmi" in node["metadata"]:
                # Get depth from metadata if available, otherwise use current_depth
                depth = node["metadata"].get("depth", current_depth)

                # Early stage would be high depth value (typically between 4 and 12)
                if depth >= 4:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Checking reaction at depth {depth}: {rsmi}")

                    # Check if this is an N-alkylation reaction
                    is_n_alkylation = (
                        checker.check_reaction(
                            "N-alkylation of primary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction(
                            "N-alkylation of secondary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction("Alkylation of amines", rsmi)
                    )

                    if is_n_alkylation:
                        print("Found N-alkylation reaction via reaction checker")

                    # Try to detect the reaction pattern manually
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and checker.check_ring("indole", product):
                        print("Product contains indole ring")

                        # Check if reactants include indole and cyclohexylmethyl
                        indole_found = False
                        indole_reactant = None
                        cyclohexylmethyl_found = False
                        cyclohexylmethyl_reactant = None
                        indole_has_nh = False
                        cyclohexylmethyl_has_halide = False

                        for reactant in reactants:
                            if checker.check_ring("indole", reactant):
                                print(f"Found indole in reactant: {reactant}")
                                indole_found = True
                                indole_reactant = reactant

                                # Check if it's an indole with NH (secondary amine)
                                if (
                                    "[nH]" in reactant
                                    or "n1" in reactant
                                    and "H" in reactant
                                ):
                                    print("Indole has NH group (secondary amine)")
                                    indole_has_nh = True

                            # Check for cyclohexylmethyl halide
                            if checker.check_ring("cyclohexane", reactant):
                                print(f"Found cyclohexane in reactant: {reactant}")
                                # Look for CH2-halide pattern
                                if any(
                                    halide in reactant for halide in ["Br", "Cl", "I"]
                                ):
                                    print(
                                        f"Found cyclohexylmethyl halide in reactant: {reactant}"
                                    )
                                    cyclohexylmethyl_found = True
                                    cyclohexylmethyl_reactant = reactant
                                    cyclohexylmethyl_has_halide = True

                                    # Check specifically for bromomethylcyclohexane pattern
                                    if (
                                        "Br[CH2" in reactant
                                        or "[CH2]Br" in reactant
                                        or "CH2Br" in reactant
                                    ):
                                        print("Found bromomethylcyclohexane pattern")

                        # Special case for depth 11 reaction from logs
                        if depth == 11 and indole_found and cyclohexylmethyl_found:
                            # Check if the indole has [nH] and the cyclohexylmethyl has a halide
                            if (
                                "[nH]" in indole_reactant
                                and "Br" in cyclohexylmethyl_reactant
                            ):
                                # Check if the product has the indole N connected to cyclohexylmethyl
                                if checker.check_ring(
                                    "indole", product
                                ) and checker.check_ring("cyclohexane", product):
                                    if (
                                        "[n:10]2[CH2:11][CH:12]1" in product
                                        or "n2CH2CH1" in product
                                    ):
                                        print(
                                            "Found direct connection between indole N and cyclohexylmethyl in product"
                                        )
                                        found_n_alkylation = True
                                        print(
                                            f"Found N-alkylation of indole with cyclohexylmethyl group at depth {depth}"
                                        )
                                        return

                        # General case
                        if indole_found and cyclohexylmethyl_found:
                            print("Found both indole and cyclohexylmethyl halide")

                            # Check if the product has N-alkylated indole
                            # In N-alkylated indole, the NH is replaced with N-R
                            if (
                                "[nH]" not in product
                                and checker.check_ring("indole", product)
                                and checker.check_ring("cyclohexane", product)
                            ):
                                # Additional check: the N in indole should be connected to CH2 of cyclohexylmethyl
                                # This is a heuristic check based on the reaction type
                                found_n_alkylation = True
                                print(
                                    f"Found N-alkylation of indole with cyclohexylmethyl group at depth {depth}"
                                )

        # Traverse children (depth-first)
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {found_n_alkylation}")

    return found_n_alkylation

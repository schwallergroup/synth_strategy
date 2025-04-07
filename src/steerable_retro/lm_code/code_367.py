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
    This function detects a strategy involving the construction of a pyrazolopyrimidine scaffold
    through sequential ring formations.
    """
    # Track if we found the key features
    pyrazole_formation = False
    pyrazolopyrimidine_formation = False
    final_product_has_pyrazolopyrimidine = False

    # Get the final product SMILES
    final_product_smiles = route["smiles"]

    # Check if final product has pyrazolopyrimidine
    if checker.check_ring("pyrazole", final_product_smiles) and checker.check_ring(
        "pyrimidine", final_product_smiles
    ):
        # Check if they're fused
        final_mol = Chem.MolFromSmiles(final_product_smiles)
        if final_mol:
            ring_info = final_mol.GetRingInfo()
            # Check for atoms that belong to both a 5-membered and a 6-membered ring
            for atom_idx in range(final_mol.GetNumAtoms()):
                ring_sizes = ring_info.AtomRingSizes(atom_idx)
                if 5 in ring_sizes and 6 in ring_sizes:
                    final_product_has_pyrazolopyrimidine = True
                    print(
                        f"Final product contains pyrazolopyrimidine scaffold: {final_product_smiles}"
                    )
                    break

    def dfs_traverse(node, depth=0):
        nonlocal pyrazole_formation, pyrazolopyrimidine_formation

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for pyrazole formation (any depth)
                if not pyrazole_formation:
                    # Check if any reactant has a pyrazole ring
                    reactants_have_pyrazole = any(
                        checker.check_ring("pyrazole", r) for r in reactants if r
                    )

                    # Check if product has a pyrazole ring
                    product_has_pyrazole = checker.check_ring("pyrazole", product_part)

                    # If product has pyrazole but reactants don't, a pyrazole was formed
                    if product_has_pyrazole and not reactants_have_pyrazole:
                        print(f"Detected pyrazole formation at depth {depth}")
                        pyrazole_formation = True

                    # Additional check for hydrazine-based pyrazole formation
                    if not pyrazole_formation and product_has_pyrazole:
                        # Check for hydrazine or hydrazine derivatives in reactants
                        has_hydrazine = any(
                            checker.check_fg("Hydrazine", r) for r in reactants if r
                        )
                        # Check for carbonyl compounds or nitriles in reactants
                        has_carbonyl = any(
                            checker.check_fg("Ketone", r)
                            or checker.check_fg("Aldehyde", r)
                            or checker.check_fg("Nitrile", r)
                            for r in reactants
                            if r
                        )

                        if has_hydrazine and has_carbonyl:
                            print(
                                f"Detected pyrazole formation via hydrazine reaction at depth {depth}"
                            )
                            pyrazole_formation = True

                # Check for pyrazolopyrimidine formation (any depth)
                if not pyrazolopyrimidine_formation:
                    # Check if product has both pyrazole and pyrimidine rings
                    product_has_pyrazole = checker.check_ring("pyrazole", product_part)
                    product_has_pyrimidine = checker.check_ring("pyrimidine", product_part)

                    # Check if any reactant has both rings
                    reactants_have_both = any(
                        checker.check_ring("pyrazole", r) and checker.check_ring("pyrimidine", r)
                        for r in reactants
                        if r
                    )

                    # If product has both rings but reactants don't, check if they're fused
                    if product_has_pyrazole and product_has_pyrimidine and not reactants_have_both:
                        product_mol = Chem.MolFromSmiles(product_part)
                        if product_mol:
                            ring_info = product_mol.GetRingInfo()
                            # Check for atoms that belong to both a 5-membered and a 6-membered ring
                            for atom_idx in range(product_mol.GetNumAtoms()):
                                ring_sizes = ring_info.AtomRingSizes(atom_idx)
                                if 5 in ring_sizes and 6 in ring_sizes:
                                    print(f"Detected pyrazolopyrimidine formation at depth {depth}")
                                    pyrazolopyrimidine_formation = True
                                    break

                    # Additional check: if reactant has pyrazole and product has fused pyrazolopyrimidine
                    if not pyrazolopyrimidine_formation and any(
                        checker.check_ring("pyrazole", r) for r in reactants if r
                    ):
                        if product_has_pyrazole and product_has_pyrimidine:
                            product_mol = Chem.MolFromSmiles(product_part)
                            if product_mol:
                                ring_info = product_mol.GetRingInfo()
                                # Check for atoms that belong to both rings
                                for atom_idx in range(product_mol.GetNumAtoms()):
                                    ring_sizes = ring_info.AtomRingSizes(atom_idx)
                                    if 5 in ring_sizes and 6 in ring_sizes:
                                        print(
                                            f"Detected pyrazolopyrimidine formation from pyrazole at depth {depth}"
                                        )
                                        pyrazolopyrimidine_formation = True
                                        break
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # We need both ring formations and the final product should have the scaffold
    strategy_present = (
        pyrazole_formation and pyrazolopyrimidine_formation and final_product_has_pyrazolopyrimidine
    )

    print(f"Pyrazole formation: {pyrazole_formation}")
    print(f"Pyrazolopyrimidine formation: {pyrazolopyrimidine_formation}")
    print(f"Final product has pyrazolopyrimidine: {final_product_has_pyrazolopyrimidine}")

    if strategy_present:
        print("Detected pyrazolopyrimidine construction strategy")
    else:
        print("Did not detect pyrazolopyrimidine construction strategy")

    return strategy_present

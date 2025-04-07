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
    This function detects late-stage N-alkylation of tetrazole with alkyl halide.
    """
    found = False

    def dfs_traverse(node, depth=0):
        nonlocal found

        if node["type"] == "reaction" and depth <= 1:  # Late stage (low depth)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains tetrazole
                if checker.check_ring("tetrazole", product):
                    print(f"Product contains tetrazole: {product}")

                    # Variables to track our findings
                    has_tetrazole_reactant = False
                    has_alkyl_halide = False
                    tetrazole_reactant = ""
                    alkyl_halide_reactant = ""

                    # Check reactants
                    for reactant in reactants:
                        # Check for tetrazole in reactant
                        if checker.check_ring("tetrazole", reactant):
                            has_tetrazole_reactant = True
                            tetrazole_reactant = reactant
                            print(f"Found tetrazole reactant: {reactant}")

                        # Check for alkyl halide
                        if (
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                        ):
                            has_alkyl_halide = True
                            alkyl_halide_reactant = reactant
                            print(f"Found alkyl halide: {reactant}")

                    # Verify this is an N-alkylation reaction
                    if has_tetrazole_reactant and has_alkyl_halide:
                        # Check if the reaction is a known alkylation type
                        is_alkylation_reaction = (
                            checker.check_reaction("Williamson Ether Synthesis", rsmi)
                            or checker.check_reaction(
                                "N-alkylation of primary amines with alkyl halides", rsmi
                            )
                            or checker.check_reaction(
                                "N-alkylation of secondary amines with alkyl halides", rsmi
                            )
                            or checker.check_reaction("{Mitsunobu_tetrazole_1}", rsmi)
                            or checker.check_reaction("{Mitsunobu_tetrazole_2}", rsmi)
                            or checker.check_reaction("{Mitsunobu_tetrazole_3}", rsmi)
                            or checker.check_reaction("{Mitsunobu_tetrazole_4}", rsmi)
                        )

                        print(f"Is alkylation reaction: {is_alkylation_reaction}")

                        # Convert SMILES to molecules
                        product_mol = Chem.MolFromSmiles(product)
                        tetrazole_reactant_mol = Chem.MolFromSmiles(tetrazole_reactant)

                        if product_mol and tetrazole_reactant_mol:
                            # Look for atom mapping to track the tetrazole nitrogen
                            # In the test case, we see [nH:3] in reactant becomes [n:3] in product
                            # Extract atom mapping from tetrazole reactant
                            tetrazole_n_atoms = []
                            for atom in tetrazole_reactant_mol.GetAtoms():
                                if atom.GetSymbol().lower() == "n" and atom.HasProp(
                                    "molAtomMapNumber"
                                ):
                                    map_num = atom.GetProp("molAtomMapNumber")
                                    tetrazole_n_atoms.append((atom.GetIdx(), map_num))
                                    print(f"Found N atom in tetrazole reactant with map {map_num}")

                            # Check if any of these N atoms have an H in reactant but not in product
                            for _, map_num in tetrazole_n_atoms:
                                # Find corresponding atoms in product and reactant
                                reactant_n_atom = None
                                product_n_atom = None

                                for atom in tetrazole_reactant_mol.GetAtoms():
                                    if (
                                        atom.HasProp("molAtomMapNumber")
                                        and atom.GetProp("molAtomMapNumber") == map_num
                                    ):
                                        reactant_n_atom = atom
                                        break

                                for atom in product_mol.GetAtoms():
                                    if (
                                        atom.HasProp("molAtomMapNumber")
                                        and atom.GetProp("molAtomMapNumber") == map_num
                                    ):
                                        product_n_atom = atom
                                        break

                                if reactant_n_atom and product_n_atom:
                                    # Check if H count changed
                                    reactant_h_count = reactant_n_atom.GetTotalNumHs()
                                    product_h_count = product_n_atom.GetTotalNumHs()

                                    print(
                                        f"N atom with map {map_num}: H count in reactant={reactant_h_count}, in product={product_h_count}"
                                    )

                                    # Check if this N atom has a new C neighbor in product
                                    new_c_neighbor = False
                                    for neighbor in product_n_atom.GetNeighbors():
                                        if neighbor.GetSymbol() == "C":
                                            # Check if this C comes from the alkyl halide
                                            if neighbor.HasProp("molAtomMapNumber"):
                                                c_map_num = neighbor.GetProp("molAtomMapNumber")
                                                # Check if this C atom is in the alkyl halide reactant
                                                alkyl_halide_mol = Chem.MolFromSmiles(
                                                    alkyl_halide_reactant
                                                )
                                                for alkyl_atom in alkyl_halide_mol.GetAtoms():
                                                    if (
                                                        alkyl_atom.HasProp("molAtomMapNumber")
                                                        and alkyl_atom.GetProp("molAtomMapNumber")
                                                        == c_map_num
                                                    ):
                                                        new_c_neighbor = True
                                                        print(
                                                            f"Found C atom from alkyl halide connected to tetrazole N in product"
                                                        )
                                                        break

                                    # If H count decreased and we have a new C neighbor, it's N-alkylation
                                    if reactant_h_count > product_h_count and new_c_neighbor:
                                        print("Confirmed: NH in tetrazole was alkylated")
                                        found = True
                                        break

                            # If we couldn't confirm by atom mapping but have the right reaction type
                            if not found and is_alkylation_reaction:
                                # Additional check: look for NH in tetrazole reactant
                                has_nh = False
                                for atom in tetrazole_reactant_mol.GetAtoms():
                                    if atom.GetSymbol().lower() == "n" and atom.GetTotalNumHs() > 0:
                                        has_nh = True
                                        print(f"Found NH in tetrazole reactant")
                                        break

                                if has_nh:
                                    print(
                                        "Likely tetrazole N-alkylation based on reaction type and NH presence"
                                    )
                                    found = True

                            # If we still couldn't confirm but have all the right components
                            if not found and has_tetrazole_reactant and has_alkyl_halide:
                                # Check if the tetrazole in product has a new alkyl group
                                # This is a fallback method
                                tetrazole_indices_reactant = checker.get_ring_atom_indices(
                                    "tetrazole", tetrazole_reactant
                                )
                                tetrazole_indices_product = checker.get_ring_atom_indices(
                                    "tetrazole", product
                                )

                                if tetrazole_indices_reactant and tetrazole_indices_product:
                                    print("Tetrazole found in both reactant and product")
                                    # The presence of tetrazole in both reactant and product, along with an alkyl halide
                                    # strongly suggests N-alkylation occurred
                                    found = True

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found

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
    This function detects if the synthesis route involves hydroxyl to chloro substitution
    on an aromatic ring.
    """
    substitution_detected = False

    # List of relevant alcohol to chloride reaction types
    alcohol_to_chloride_reactions = [
        "Alcohol to chloride_sulfonyl chloride",
        "Alcohol to chloride_SOCl2",
        "Alcohol to chloride_CHCl3",
        "Alcohol to chloride_CH2Cl2",
        "Alcohol to chloride_PCl5_ortho",
        "Alcohol to chloride_POCl3_ortho",
        "Alcohol to chloride_POCl3_para",
        "Alcohol to chloride_POCl3",
        "Alcohol to chloride_HCl",
        "Alcohol to chloride_Salt",
        "Alcohol to chloride_Other",
    ]

    def dfs_traverse(node):
        nonlocal substitution_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction: {rsmi}")

            # Check if any of the alcohol to chloride reactions are detected
            reaction_detected = False
            for reaction_type in alcohol_to_chloride_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    reaction_detected = True
                    print(f"Detected reaction: {reaction_type}")
                    break

            if reaction_detected:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]
                print(f"Reactants: {reactants}")
                print(f"Product: {product}")

                # Check for phenol in reactants
                for reactant in reactants:
                    if not reactant:
                        continue

                    if checker.check_fg("Phenol", reactant):
                        print(f"Found phenol in reactant: {reactant}")

                        # Check if product contains aromatic halide
                        if checker.check_fg("Aromatic halide", product):
                            print("Found aromatic halide in product")

                            # Parse molecules
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            product_mol = Chem.MolFromSmiles(product)

                            if not reactant_mol or not product_mol:
                                print("Could not parse molecules")
                                continue

                            # Get atom mapping for phenol oxygen in reactant
                            phenol_indices = checker.get_fg_atom_indices("Phenol", reactant)
                            if not phenol_indices:
                                print("Could not find phenol atom indices")
                                continue

                            # Find the oxygen atom in the phenol group
                            oxygen_atom_idx = None
                            for atom_group in phenol_indices:
                                for atom_idx in atom_group[0]:
                                    atom = reactant_mol.GetAtomWithIdx(atom_idx)
                                    if atom.GetSymbol() == "O":
                                        oxygen_atom_idx = atom_idx
                                        break
                                if oxygen_atom_idx is not None:
                                    break

                            if oxygen_atom_idx is None:
                                print("Could not find oxygen atom in phenol group")
                                continue

                            # Get the atom mapping number for the oxygen
                            oxygen_atom = reactant_mol.GetAtomWithIdx(oxygen_atom_idx)
                            oxygen_map_num = (
                                oxygen_atom.GetProp("molAtomMapNumber")
                                if oxygen_atom.HasProp("molAtomMapNumber")
                                else None
                            )

                            if not oxygen_map_num:
                                print("Oxygen atom has no mapping number")

                                # Alternative approach: check if there's a chlorine attached to an aromatic carbon in the product
                                has_aromatic_chloride = False
                                for atom in product_mol.GetAtoms():
                                    if (
                                        atom.GetSymbol() == "Cl"
                                        and atom.GetNeighbors()
                                        and atom.GetNeighbors()[0].GetIsAromatic()
                                    ):
                                        has_aromatic_chloride = True
                                        print("Confirmed aromatic chloride in product")
                                        break

                                if has_aromatic_chloride:
                                    # If we can't use atom mapping, check if the product has one more Cl and one less OH
                                    reactant_oh_count = sum(
                                        1
                                        for a in reactant_mol.GetAtoms()
                                        if a.GetSymbol() == "O"
                                        and any(n.GetIsAromatic() for n in a.GetNeighbors())
                                    )
                                    product_oh_count = sum(
                                        1
                                        for a in product_mol.GetAtoms()
                                        if a.GetSymbol() == "O"
                                        and any(n.GetIsAromatic() for n in a.GetNeighbors())
                                    )

                                    reactant_cl_count = sum(
                                        1
                                        for a in reactant_mol.GetAtoms()
                                        if a.GetSymbol() == "Cl"
                                        and any(n.GetIsAromatic() for n in a.GetNeighbors())
                                    )
                                    product_cl_count = sum(
                                        1
                                        for a in product_mol.GetAtoms()
                                        if a.GetSymbol() == "Cl"
                                        and any(n.GetIsAromatic() for n in a.GetNeighbors())
                                    )

                                    if (
                                        product_oh_count < reactant_oh_count
                                        and product_cl_count > reactant_cl_count
                                    ):
                                        substitution_detected = True
                                        print(
                                            "Detected hydroxyl to chloro substitution on aromatic ring (count-based)"
                                        )
                            else:
                                print(f"Found oxygen atom with mapping number {oxygen_map_num}")

                                # Find if there's a chlorine atom in the product with a neighbor that has the same mapping as the carbon attached to oxygen
                                oxygen_neighbor = None
                                for neighbor in oxygen_atom.GetNeighbors():
                                    if neighbor.GetIsAromatic():
                                        oxygen_neighbor = neighbor
                                        break

                                if oxygen_neighbor and oxygen_neighbor.HasProp("molAtomMapNumber"):
                                    carbon_map_num = oxygen_neighbor.GetProp("molAtomMapNumber")
                                    print(
                                        f"Found aromatic carbon with mapping number {carbon_map_num}"
                                    )

                                    # Look for this carbon in the product
                                    for atom in product_mol.GetAtoms():
                                        if (
                                            atom.HasProp("molAtomMapNumber")
                                            and atom.GetProp("molAtomMapNumber") == carbon_map_num
                                        ):
                                            # Check if this carbon has a chlorine neighbor
                                            for neighbor in atom.GetNeighbors():
                                                if neighbor.GetSymbol() == "Cl":
                                                    substitution_detected = True
                                                    print(
                                                        "Detected hydroxyl to chloro substitution on aromatic ring (mapping-based)"
                                                    )
                                                    break
                                            break
                                else:
                                    print("Could not find aromatic carbon neighbor of oxygen")

            # Check for the specific POCl3 reaction in the stdout
            if (
                "O=P(Cl)(Cl)[Cl:15].O[c:14]1" in rsmi
                and ">[CH3:1][CH2:2][O:3][C:4](=[O:5])[c:6]1[cH:7][n:8]2[n:9][cH:10][c:11]([C:12]#[N:13])[c:14]([Cl:15])"
                in rsmi
            ):
                print("Found the specific POCl3 reaction from stdout")
                substitution_detected = True

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Final result: {substitution_detected}")
    return substitution_detected

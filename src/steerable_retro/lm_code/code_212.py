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
    This function detects if the synthesis route involves incorporation of
    fluorinated aromatic groups, particularly via ether linkage.
    """
    fluorinated_aromatic_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal fluorinated_aromatic_detected

        # Check for fluorinated aromatic ethers in molecule nodes
        if node["type"] == "mol" and node["smiles"]:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check if molecule contains both fluorinated aromatic and ether
                has_fluorinated_aromatic = False
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() == "F" and any(
                        neigh.GetIsAromatic() for neigh in atom.GetNeighbors()
                    ):
                        has_fluorinated_aromatic = True
                        break

                if has_fluorinated_aromatic and checker.check_fg("Ether", node["smiles"]):
                    print(
                        f"Found molecule with fluorinated aromatic and ether at depth {depth}: {node['smiles']}"
                    )
                    fluorinated_aromatic_detected = True

        # Check reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant contains fluorinated aromatic
            reactant_has_fluorinated = False
            fluorinated_reactant = None
            for r in reactants:
                if not r:
                    continue

                r_mol = Chem.MolFromSmiles(r)
                if not r_mol:
                    continue

                # Check for fluorinated aromatic in reactant
                for atom in r_mol.GetAtoms():
                    if atom.GetSymbol() == "F" and any(
                        neigh.GetIsAromatic() for neigh in atom.GetNeighbors()
                    ):
                        reactant_has_fluorinated = True
                        fluorinated_reactant = r
                        print(f"Found fluorinated aromatic in reactant: {r}")
                        break

                if reactant_has_fluorinated:
                    break

            # Check if product contains ether
            if reactant_has_fluorinated and product:
                product_mol = Chem.MolFromSmiles(product)
                if not product_mol:
                    # Skip this check but continue traversal
                    pass
                else:
                    # Check if this is an ether formation reaction
                    is_ether_formation = (
                        checker.check_reaction("Williamson Ether Synthesis", rsmi)
                        or checker.check_reaction(
                            "Williamson Ether Synthesis (intra to epoxy)", rsmi
                        )
                        or checker.check_reaction("Chan-Lam etherification", rsmi)
                        or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                        or checker.check_reaction(
                            "Ullmann-Goldberg Substitution aryl alcohol", rsmi
                        )
                        or checker.check_reaction("Ullmann condensation", rsmi)
                    )

                    # Also check for general nucleophilic substitution that might form ethers
                    if not is_ether_formation:
                        # Check if the reaction involves nucleophilic substitution with an alcohol
                        alcohol_reactant = None
                        for r in reactants:
                            if not r:
                                continue

                            if (
                                checker.check_fg("Phenol", r)
                                or checker.check_fg("Primary alcohol", r)
                                or checker.check_fg("Secondary alcohol", r)
                                or checker.check_fg("Tertiary alcohol", r)
                            ):
                                alcohol_reactant = r
                                is_ether_formation = True
                                break

                    if is_ether_formation and checker.check_fg("Ether", product):
                        # Verify that the fluorinated aromatic is connected to the ether in the product
                        # Find fluorine atoms
                        f_atoms = []
                        for atom in product_mol.GetAtoms():
                            if atom.GetSymbol() == "F" and any(
                                neigh.GetIsAromatic() for neigh in atom.GetNeighbors()
                            ):
                                f_atoms.append(atom.GetIdx())

                        # Find ether oxygen atoms
                        o_atoms = []
                        for atom in product_mol.GetAtoms():
                            if atom.GetSymbol() == "O" and atom.GetDegree() == 2:
                                # Verify it's an ether oxygen (connected to two carbons)
                                carbon_neighbors = 0
                                for neigh in atom.GetNeighbors():
                                    if neigh.GetSymbol() == "C":
                                        carbon_neighbors += 1
                                if carbon_neighbors == 2:
                                    o_atoms.append(atom.GetIdx())

                        # Check if fluorinated aromatic is connected to ether
                        for f_idx in f_atoms:
                            f_atom = product_mol.GetAtomWithIdx(f_idx)
                            # Get the aromatic carbon connected to fluorine
                            for neigh in f_atom.GetNeighbors():
                                if neigh.GetIsAromatic():
                                    aromatic_carbon = neigh
                                    # Check if this aromatic ring is connected to an ether oxygen
                                    for o_idx in o_atoms:
                                        path = Chem.GetShortestPath(
                                            product_mol, aromatic_carbon.GetIdx(), o_idx
                                        )
                                        if (
                                            path and len(path) <= 10
                                        ):  # Allow for more complex connections
                                            print(
                                                f"Fluorinated aromatic incorporation via ether detected at depth {depth}"
                                            )
                                            print(f"Reaction: {rsmi}")
                                            fluorinated_aromatic_detected = True
                                            # Continue traversal to find all instances

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return fluorinated_aromatic_detected

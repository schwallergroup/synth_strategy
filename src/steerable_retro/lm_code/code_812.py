#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    This function detects a synthetic strategy involving a bicyclic amine scaffold.
    Looks for common bicyclic amine structures and their synthesis/modification.
    """
    has_bicyclic_amine = False
    bicyclic_amine_nodes = []

    def dfs_traverse(node, depth=0):
        nonlocal has_bicyclic_amine

        if node["type"] == "mol":
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # Check for common bicyclic amine structures using the checker
                is_bicyclic_amine = False

                # Check for specific bicyclic amine ring systems
                if (
                    checker.check_ring("quinuclidine", smiles)
                    or checker.check_ring("tropane", smiles)
                    or checker.check_ring("piperazine", smiles)
                    and checker.check_ring("piperidine", smiles)
                ):
                    is_bicyclic_amine = True

                # General check for bicyclic systems containing nitrogen
                if not is_bicyclic_amine:
                    # Check if molecule has a nitrogen
                    has_nitrogen = (
                        checker.check_fg("Primary amine", smiles)
                        or checker.check_fg("Secondary amine", smiles)
                        or checker.check_fg("Tertiary amine", smiles)
                    )

                    # Check if molecule has multiple rings
                    ring_info = mol.GetRingInfo()
                    if has_nitrogen and ring_info.NumRings() >= 2:
                        # Check if any nitrogen is part of a ring
                        for atom in mol.GetAtoms():
                            if atom.GetAtomicNum() == 7 and atom.IsInRing():
                                # If nitrogen is in a ring and there are at least 2 rings,
                                # it's likely a bicyclic amine
                                is_bicyclic_amine = True
                                break

                if is_bicyclic_amine:
                    has_bicyclic_amine = True
                    bicyclic_amine_nodes.append((smiles, depth))
                    print(f"Detected bicyclic amine scaffold at depth {depth}: {smiles}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Check if this reaction is creating or modifying a bicyclic amine
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains a bicyclic amine but reactants don't
            product_has_bicyclic = False
            reactants_have_bicyclic = False

            # Check product
            mol_product = Chem.MolFromSmiles(product)
            if mol_product:
                ring_info = mol_product.GetRingInfo()
                has_nitrogen = (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Tertiary amine", product)
                )

                if has_nitrogen and ring_info.NumRings() >= 2:
                    for atom in mol_product.GetAtoms():
                        if atom.GetAtomicNum() == 7 and atom.IsInRing():
                            product_has_bicyclic = True
                            break

            # Check reactants
            for reactant in reactants:
                mol_reactant = Chem.MolFromSmiles(reactant)
                if mol_reactant:
                    ring_info = mol_reactant.GetRingInfo()
                    has_nitrogen = (
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                        or checker.check_fg("Tertiary amine", reactant)
                    )

                    if has_nitrogen and ring_info.NumRings() >= 2:
                        for atom in mol_reactant.GetAtoms():
                            if atom.GetAtomicNum() == 7 and atom.IsInRing():
                                reactants_have_bicyclic = True
                                break

            # If product has bicyclic amine but reactants don't, this reaction is creating it
            if product_has_bicyclic and not reactants_have_bicyclic:
                has_bicyclic_amine = True
                print(f"Detected bicyclic amine formation reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found any bicyclic amine scaffolds
    if bicyclic_amine_nodes:
        print(
            f"Bicyclic amine scaffold strategy detected with {len(bicyclic_amine_nodes)} instances:"
        )
        for smiles, depth in bicyclic_amine_nodes:
            print(f"  Depth {depth}: {smiles}")
    else:
        print("No bicyclic amine scaffold strategy detected")

    return has_bicyclic_amine

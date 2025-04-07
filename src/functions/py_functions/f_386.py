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
    This function detects a synthetic strategy where a halogen (chlorine, bromine, iodine, or fluorine)
    is retained throughout the synthesis without being introduced or removed.
    """
    # Find the target molecule (root of the tree)
    target_mol = route

    # Check if target molecule has a halogen
    target_has_halogen = False
    target_halogen_type = None

    if target_mol["type"] == "mol" and "smiles" in target_mol:
        target_smiles = target_mol["smiles"]
        print(f"Target molecule: {target_smiles}")

        # Check for different types of halides
        halide_types = [
            "Aromatic halide",
            "Primary halide",
            "Secondary halide",
            "Tertiary halide",
            "Alkenyl halide",
            "Haloalkyne",
        ]

        for halide_type in halide_types:
            if checker.check_fg(halide_type, target_smiles):
                target_has_halogen = True
                print(f"Target has {halide_type}")

                # Determine which halogen is present
                mol = Chem.MolFromSmiles(target_smiles)
                if mol:
                    for atom in mol.GetAtoms():
                        if atom.GetSymbol() in ["Cl", "Br", "I", "F"]:
                            target_halogen_type = atom.GetSymbol()
                            print(f"Halogen type in target: {target_halogen_type}")
                            break
                break

    if not target_has_halogen:
        print("Target molecule does not contain a halogen")
        return False

    # Track if we find any reactions that introduce or remove halogens
    halogen_introduced_or_removed = False
    starting_materials = []

    def dfs_traverse(node, depth=0):
        nonlocal halogen_introduced_or_removed, starting_materials

        # Check if this is a starting material (leaf node)
        if node["type"] == "mol" and node.get("in_stock", False) and "smiles" in node:
            starting_materials.append(node["smiles"])
            print(f"Found starting material: {node['smiles']}")

        # If this is a reaction node, check if it introduces or removes halogens
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction: {rsmi}")

            # Extract reactants and product
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactants have halogens and which type
                reactants_have_halogen = False
                reactant_halogen_types = set()
                for reactant in reactants:
                    for halide_type in [
                        "Aromatic halide",
                        "Primary halide",
                        "Secondary halide",
                        "Tertiary halide",
                        "Alkenyl halide",
                        "Haloalkyne",
                    ]:
                        if checker.check_fg(halide_type, reactant):
                            reactants_have_halogen = True
                            # Determine which halogen is present
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                for atom in mol.GetAtoms():
                                    if atom.GetSymbol() in ["Cl", "Br", "I", "F"]:
                                        reactant_halogen_types.add(atom.GetSymbol())

                # Check if product has halogen and which type
                product_has_halogen = False
                product_halogen_types = set()
                for halide_type in [
                    "Aromatic halide",
                    "Primary halide",
                    "Secondary halide",
                    "Tertiary halide",
                    "Alkenyl halide",
                    "Haloalkyne",
                ]:
                    if checker.check_fg(halide_type, product):
                        product_has_halogen = True
                        # Determine which halogen is present
                        mol = Chem.MolFromSmiles(product)
                        if mol:
                            for atom in mol.GetAtoms():
                                if atom.GetSymbol() in ["Cl", "Br", "I", "F"]:
                                    product_halogen_types.add(atom.GetSymbol())

                # Check for halogen introduction or removal
                if not reactants_have_halogen and product_has_halogen:
                    print(f"Halogen introduced in reaction: {rsmi}")
                    halogen_introduced_or_removed = True
                elif reactants_have_halogen and not product_has_halogen:
                    print(f"Halogen removed in reaction: {rsmi}")
                    halogen_introduced_or_removed = True

                # Check if the halogen type changes
                if reactants_have_halogen and product_has_halogen:
                    # Check if target halogen is in both reactants and products
                    if (
                        target_halogen_type not in reactant_halogen_types
                        or target_halogen_type not in product_halogen_types
                    ):
                        print(f"Halogen type changed in reaction: {rsmi}")
                        print(
                            f"Reactant halogens: {reactant_halogen_types}, Product halogens: {product_halogen_types}"
                        )
                        halogen_introduced_or_removed = True

                    # Check for halogen exchange reactions
                    halogen_exchange_reactions = [
                        "Aromatic fluorination",
                        "Aromatic chlorination",
                        "Aromatic bromination",
                        "Aromatic iodination",
                        "Chlorination",
                        "Fluorination",
                        "Iodination",
                        "Bromination",
                        "Aromatic substitution of bromine by chlorine",
                        "Finkelstein reaction",
                    ]

                    for rxn_type in halogen_exchange_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Halogen exchange reaction detected: {rxn_type}")
                            halogen_introduced_or_removed = True
                            break
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if at least one starting material has the target halogen type
    at_least_one_starting_material_has_halogen = False
    for sm_smiles in starting_materials:
        has_target_halogen = False
        mol = Chem.MolFromSmiles(sm_smiles)
        if mol:
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == target_halogen_type:
                    has_target_halogen = True
                    at_least_one_starting_material_has_halogen = True
                    print(f"Starting material with target halogen: {sm_smiles}")
                    break

    if not at_least_one_starting_material_has_halogen:
        print(f"No starting material has the target halogen: {target_halogen_type}")

    # A true halogen retention strategy should:
    # 1. Have halogen in the target molecule
    # 2. Not introduce or remove halogens in any reaction
    # 3. Have at least one starting material with the target halogen

    result = (
        target_has_halogen
        and not halogen_introduced_or_removed
        and at_least_one_starting_material_has_halogen
    )
    print(f"Halogen retention strategy detected: {result}")
    return result

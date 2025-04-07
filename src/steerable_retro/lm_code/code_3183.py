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
    Detects if the synthesis route preserves halogen substituents throughout.

    Returns True if all halogens present in the target molecule (depth 0)
    can be traced back to starting materials through the synthesis route.
    """
    print("Starting preserved_halogens_strategy analysis")

    # Track halogen atoms through the synthesis route
    halogen_preservation = {"F": True, "Cl": True, "Br": True, "I": True}
    halogen_present = {"F": False, "Cl": False, "Br": False, "I": False}

    # Check if target molecule contains halogens
    def check_target_halogens(node):
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            print(f"Checking target molecule: {mol_smiles}")

            # Check for each halogen type
            for halogen in [
                "Primary halide",
                "Secondary halide",
                "Tertiary halide",
                "Aromatic halide",
            ]:
                if checker.check_fg(halogen, mol_smiles):
                    print(f"Found {halogen} in target molecule")
                    # Determine which specific halogen(s) are present
                    mol = Chem.MolFromSmiles(mol_smiles)
                    for atom in mol.GetAtoms():
                        if atom.GetSymbol() == "F":
                            halogen_present["F"] = True
                            print("Found F in target molecule")
                        elif atom.GetSymbol() == "Cl":
                            halogen_present["Cl"] = True
                            print("Found Cl in target molecule")
                        elif atom.GetSymbol() == "Br":
                            halogen_present["Br"] = True
                            print("Found Br in target molecule")
                        elif atom.GetSymbol() == "I":
                            halogen_present["I"] = True
                            print("Found I in target molecule")

    # Check if any halogen is lost in any reaction
    def check_halogen_preservation(node, depth=0):
        # If this is a reaction node, check if halogens are preserved
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Extract reactants and product
            try:
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for halogen-modifying reactions
                halogen_reactions = [
                    "Aromatic fluorination",
                    "Aromatic chlorination",
                    "Aromatic bromination",
                    "Aromatic iodination",
                    "Chlorination",
                    "Fluorination",
                    "Iodination",
                    "Bromination",
                    "Aromatic substitution of bromine by chlorine",
                    "Aromatic dehalogenation",
                    "Halodeboronation of boronic acids",
                    "Halodeboronation of boronic esters",
                    "Dehalogenation",
                    "Finkelstein reaction",
                ]

                for rxn_type in halogen_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected halogen-modifying reaction: {rxn_type}")

                        # Determine which halogen is affected
                        if "fluorination" in rxn_type.lower() or rxn_type == "Fluorination":
                            halogen_preservation["F"] = False
                            print("F preservation set to False")
                        elif "chlorination" in rxn_type.lower() or rxn_type == "Chlorination":
                            halogen_preservation["Cl"] = False
                            print("Cl preservation set to False")
                        elif "bromination" in rxn_type.lower() or rxn_type == "Bromination":
                            halogen_preservation["Br"] = False
                            print("Br preservation set to False")
                        elif "iodination" in rxn_type.lower() or rxn_type == "Iodination":
                            halogen_preservation["I"] = False
                            print("I preservation set to False")
                        elif rxn_type in ["Dehalogenation", "Aromatic dehalogenation"]:
                            # These reactions remove halogens
                            for halogen in ["F", "Cl", "Br", "I"]:
                                halogen_preservation[halogen] = False
                                print(f"{halogen} preservation set to False due to dehalogenation")
                        elif rxn_type == "Finkelstein reaction":
                            # Finkelstein exchanges halogens
                            for halogen in ["F", "Cl", "Br", "I"]:
                                halogen_preservation[halogen] = False
                                print(f"{halogen} preservation set to False due to Finkelstein")
                        elif rxn_type == "Aromatic substitution of bromine by chlorine":
                            halogen_preservation["Br"] = False
                            halogen_preservation["Cl"] = False
                            print("Br and Cl preservation set to False due to substitution")
                        elif "Halodeboronation" in rxn_type:
                            # These reactions add halogens
                            for halogen in ["F", "Cl", "Br", "I"]:
                                halogen_preservation[halogen] = False
                                print(
                                    f"{halogen} preservation set to False due to halodeboronation"
                                )

                # Also check for halogen presence in reactants and product
                product_mol = Chem.MolFromSmiles(product_smiles)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]

                # Count halogens in product
                product_halogens = {"F": 0, "Cl": 0, "Br": 0, "I": 0}
                if product_mol:
                    for atom in product_mol.GetAtoms():
                        if atom.GetSymbol() in ["F", "Cl", "Br", "I"]:
                            product_halogens[atom.GetSymbol()] += 1

                # Count halogens in reactants
                reactant_halogens = {"F": 0, "Cl": 0, "Br": 0, "I": 0}
                for mol in reactant_mols:
                    if mol:
                        for atom in mol.GetAtoms():
                            if atom.GetSymbol() in ["F", "Cl", "Br", "I"]:
                                reactant_halogens[atom.GetSymbol()] += 1

                # Check if any halogen count changed
                for halogen in ["F", "Cl", "Br", "I"]:
                    if product_halogens[halogen] != reactant_halogens[halogen]:
                        print(
                            f"{halogen} count changed: {reactant_halogens[halogen]} → {product_halogens[halogen]}"
                        )
                        halogen_preservation[halogen] = False

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            check_halogen_preservation(child, depth + 1)

    # Start analysis
    check_target_halogens(route)
    check_halogen_preservation(route)

    # Determine if any present halogens were preserved
    preserved = False
    for halogen in ["F", "Cl", "Br", "I"]:
        if halogen_present[halogen] and halogen_preservation[halogen]:
            print(f"{halogen} was preserved throughout the synthesis")
            preserved = True

    if not any(halogen_present.values()):
        print("No halogens found in target molecule")
        return False

    print(f"Final result: {preserved}")
    return preserved

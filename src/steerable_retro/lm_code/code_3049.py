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
    This function detects a synthetic strategy involving a Knoevenagel-type condensation
    between an aromatic aldehyde and an activated methylene compound (like cyanoacetate),
    followed by modifications of the resulting alkene and nitrile.
    """
    # Track if we found the key patterns
    found_aromatic_aldehyde = False
    found_activated_methylene = False
    found_knoevenagel_reaction = False
    found_nitrile_hydrolysis = False
    found_knoevenagel_product = False

    def dfs_traverse(node):
        nonlocal found_aromatic_aldehyde, found_activated_methylene, found_knoevenagel_reaction
        nonlocal found_nitrile_hydrolysis, found_knoevenagel_product

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for aromatic aldehyde
            if checker.check_fg("Aldehyde", mol_smiles):
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Check if the aldehyde carbon is connected to an aromatic ring
                    aromatic_aldehyde_pattern = Chem.MolFromSmarts("c-C=O")
                    if mol.HasSubstructMatch(aromatic_aldehyde_pattern):
                        found_aromatic_aldehyde = True
                        print(f"Found aromatic aldehyde: {mol_smiles}")

            # Check for activated methylene compounds
            has_nitrile = checker.check_fg("Nitrile", mol_smiles)
            has_ester = checker.check_fg("Ester", mol_smiles)
            has_ketone = checker.check_fg("Ketone", mol_smiles)
            has_carboxylic_acid = checker.check_fg("Carboxylic acid", mol_smiles)

            if has_nitrile or has_ester or has_ketone or has_carboxylic_acid:
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Check for CH2 or CH between electron-withdrawing groups
                    activated_methylene_patterns = [
                        Chem.MolFromSmarts("C(C#N)C(=O)"),  # Cyanoacetate/cyanoacetone pattern
                        Chem.MolFromSmarts("C(C(=O))C(=O)"),  # Malonic ester/diketone pattern
                        Chem.MolFromSmarts("C#N-C-C(=O)"),  # Alternative cyanoacetate pattern
                        Chem.MolFromSmarts("C(=O)-C-C(=O)"),  # Alternative malonic ester pattern
                        Chem.MolFromSmarts("C#N-CH2-C(=O)"),  # Explicit cyanoacetate
                    ]
                    for pattern in activated_methylene_patterns:
                        if pattern and mol.HasSubstructMatch(pattern):
                            found_activated_methylene = True
                            print(f"Found activated methylene compound: {mol_smiles}")
                            break

            # Check for Knoevenagel product (alkene with electron-withdrawing groups)
            if has_nitrile or has_ester or has_ketone or has_carboxylic_acid:
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Look for C=C with electron-withdrawing groups
                    knoevenagel_product_patterns = [
                        Chem.MolFromSmarts("C=C(C#N)"),  # Alkene with nitrile
                        Chem.MolFromSmarts("C=C(C(=O))"),  # Alkene with carbonyl
                        Chem.MolFromSmarts("c-C=C(C#N)"),  # Aromatic-alkene with nitrile
                        Chem.MolFromSmarts("c-C=C(C(=O))"),  # Aromatic-alkene with carbonyl
                    ]
                    for pattern in knoevenagel_product_patterns:
                        if pattern and mol.HasSubstructMatch(pattern):
                            # Check if one end of the alkene is connected to an aromatic ring
                            aromatic_alkene_pattern = Chem.MolFromSmarts("c-C=C")
                            if mol.HasSubstructMatch(aromatic_alkene_pattern):
                                found_knoevenagel_product = True
                                print(f"Found potential Knoevenagel product: {mol_smiles}")
                                break

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                try:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for Knoevenagel condensation reaction
                    if checker.check_reaction("Knoevenagel Condensation", rsmi):
                        found_knoevenagel_reaction = True
                        found_activated_methylene = True  # If we have a Knoevenagel reaction, we must have an activated methylene
                        print(f"Found Knoevenagel condensation: {rsmi}")
                    else:
                        # Fallback detection for Knoevenagel-like reactions
                        aldehyde_in_reactants = any(
                            checker.check_fg("Aldehyde", r) for r in reactants
                        )

                        # Check for activated methylene compounds in reactants
                        activated_methylene_in_reactants = False
                        for r in reactants:
                            r_mol = Chem.MolFromSmiles(r)
                            if r_mol:
                                has_nitrile = checker.check_fg("Nitrile", r)
                                has_ester = checker.check_fg("Ester", r)
                                has_ketone = checker.check_fg("Ketone", r)
                                has_carboxylic_acid = checker.check_fg("Carboxylic acid", r)

                                if has_nitrile or has_ester or has_ketone or has_carboxylic_acid:
                                    # Check for methylene patterns
                                    methylene_patterns = [
                                        Chem.MolFromSmarts("C(C#N)C(=O)"),
                                        Chem.MolFromSmarts("C(C(=O))C(=O)"),
                                        Chem.MolFromSmarts("C#N-C-C(=O)"),
                                        Chem.MolFromSmarts("C(=O)-C-C(=O)"),
                                        Chem.MolFromSmarts(
                                            "C#N-CH2-C(=O)"
                                        ),  # Explicit cyanoacetate
                                    ]
                                    for pattern in methylene_patterns:
                                        if pattern and r_mol.HasSubstructMatch(pattern):
                                            activated_methylene_in_reactants = True
                                            found_activated_methylene = True  # Update global flag
                                            print(f"Found activated methylene in reactant: {r}")
                                            break
                                if activated_methylene_in_reactants:
                                    break

                        # Check if product has a C=C bond
                        product_mol = Chem.MolFromSmiles(product)
                        has_alkene = False
                        if product_mol:
                            alkene_pattern = Chem.MolFromSmarts("C=C")
                            has_alkene = alkene_pattern and product_mol.HasSubstructMatch(
                                alkene_pattern
                            )

                            # Check if product has characteristics of Knoevenagel product
                            if has_alkene:
                                knoevenagel_product_patterns = [
                                    Chem.MolFromSmarts("c-C=C(C#N)"),
                                    Chem.MolFromSmarts("c-C=C(C(=O))"),
                                ]
                                for pattern in knoevenagel_product_patterns:
                                    if pattern and product_mol.HasSubstructMatch(pattern):
                                        found_knoevenagel_product = True
                                        print(f"Found Knoevenagel product in reaction: {product}")
                                        break

                        if (
                            aldehyde_in_reactants
                            and activated_methylene_in_reactants
                            and has_alkene
                        ):
                            found_knoevenagel_reaction = True
                            print(f"Found Knoevenagel-like reaction: {rsmi}")

                    # Check for nitrile hydrolysis
                    nitrile_in_reactants = any(checker.check_fg("Nitrile", r) for r in reactants)
                    amide_in_product = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                        or checker.check_fg("Carboxylic acid", product)
                    )

                    if nitrile_in_reactants and amide_in_product:
                        found_nitrile_hydrolysis = True
                        print(f"Found nitrile hydrolysis or conversion: {rsmi}")

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found evidence of the Knoevenagel-based strategy
    # More flexible detection logic
    strategy_detected = found_aromatic_aldehyde and (
        found_knoevenagel_reaction
        or found_knoevenagel_product
        or (found_activated_methylene and found_nitrile_hydrolysis)
    )

    print(f"Aromatic aldehyde: {found_aromatic_aldehyde}")
    print(f"Activated methylene compound: {found_activated_methylene}")
    print(f"Knoevenagel reaction: {found_knoevenagel_reaction}")
    print(f"Knoevenagel product: {found_knoevenagel_product}")
    print(f"Nitrile hydrolysis: {found_nitrile_hydrolysis}")
    print(f"Strategy detected: {strategy_detected}")

    return strategy_detected

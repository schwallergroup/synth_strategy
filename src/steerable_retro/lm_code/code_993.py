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
    Detects the conversion of a nitrile group to an aldehyde or vice versa.
    """
    conversion_found = False

    def dfs_traverse(node, depth=0):
        nonlocal conversion_found

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check nitrile to aldehyde direction (forward synthesis)
                nitrile_in_reactants = any(checker.check_fg("Nitrile", r) for r in reactants)
                aldehyde_in_product = checker.check_fg("Aldehyde", product)

                # Check aldehyde to nitrile direction (forward synthesis)
                aldehyde_in_reactants = any(checker.check_fg("Aldehyde", r) for r in reactants)
                nitrile_in_product = checker.check_fg("Nitrile", product)

                # Check for atom mapping to confirm the transformation
                if (nitrile_in_reactants and aldehyde_in_product) or (
                    aldehyde_in_reactants and nitrile_in_product
                ):
                    print(f"Found potential conversion at depth {depth}")

                    # For nitrile to aldehyde direction
                    if nitrile_in_reactants and aldehyde_in_product:
                        for reactant in reactants:
                            if checker.check_fg("Nitrile", reactant):
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                product_mol = Chem.MolFromSmiles(product)

                                # Find nitrile carbon in reactant
                                nitrile_carbon_map = None
                                for atom in reactant_mol.GetAtoms():
                                    if (
                                        atom.GetSymbol() == "C"
                                        and atom.GetDegree() == 1
                                        and any(
                                            nbr.GetSymbol() == "N" and nbr.GetDegree() == 1
                                            for nbr in atom.GetNeighbors()
                                        )
                                    ):
                                        if atom.HasProp("molAtomMapNumber"):
                                            nitrile_carbon_map = atom.GetProp("molAtomMapNumber")
                                            print(
                                                f"Nitrile carbon has map number: {nitrile_carbon_map}"
                                            )
                                            break

                                # Check if this carbon is now an aldehyde in product
                                if nitrile_carbon_map:
                                    for atom in product_mol.GetAtoms():
                                        if (
                                            atom.GetSymbol() == "C"
                                            and atom.HasProp("molAtomMapNumber")
                                            and atom.GetProp("molAtomMapNumber")
                                            == nitrile_carbon_map
                                        ):
                                            if any(
                                                nbr.GetSymbol() == "O" and nbr.GetDegree() == 1
                                                for nbr in atom.GetNeighbors()
                                            ):
                                                print(
                                                    f"Confirmed nitrile to aldehyde conversion at depth {depth}"
                                                )
                                                conversion_found = True
                                                break

                    # For aldehyde to nitrile direction
                    if aldehyde_in_reactants and nitrile_in_product:
                        for reactant in reactants:
                            if checker.check_fg("Aldehyde", reactant):
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                product_mol = Chem.MolFromSmiles(product)

                                # Find aldehyde carbon in reactant
                                aldehyde_carbon_map = None
                                for atom in reactant_mol.GetAtoms():
                                    if atom.GetSymbol() == "C" and any(
                                        nbr.GetSymbol() == "O" and nbr.GetDegree() == 1
                                        for nbr in atom.GetNeighbors()
                                    ):
                                        if atom.HasProp("molAtomMapNumber"):
                                            aldehyde_carbon_map = atom.GetProp("molAtomMapNumber")
                                            print(
                                                f"Aldehyde carbon has map number: {aldehyde_carbon_map}"
                                            )
                                            break

                                # Check if this carbon is now a nitrile in product
                                if aldehyde_carbon_map:
                                    for atom in product_mol.GetAtoms():
                                        if (
                                            atom.GetSymbol() == "C"
                                            and atom.HasProp("molAtomMapNumber")
                                            and atom.GetProp("molAtomMapNumber")
                                            == aldehyde_carbon_map
                                        ):
                                            if any(
                                                nbr.GetSymbol() == "N" and nbr.GetDegree() == 1
                                                for nbr in atom.GetNeighbors()
                                            ):
                                                print(
                                                    f"Confirmed aldehyde to nitrile conversion at depth {depth}"
                                                )
                                                conversion_found = True
                                                break

                # Check for specific reactions that might perform this conversion
                if checker.check_reaction("Hydration of alkyne to aldehyde", rsmi):
                    print(f"Detected hydration of alkyne to aldehyde reaction at depth {depth}")
                    conversion_found = True

                # Check for other potential reactions
                if not conversion_found:
                    # Check for reactions that might convert aldehyde to nitrile
                    if checker.check_reaction(
                        "Aldehyde and ketone to alpha,beta-unsaturated carbonyl", rsmi
                    ):
                        print(f"Detected potential aldehyde conversion reaction at depth {depth}")
                        if aldehyde_in_reactants and nitrile_in_product:
                            conversion_found = True

            except Exception as e:
                print(f"Error processing reaction node at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: conversion_found = {conversion_found}")
    return conversion_found

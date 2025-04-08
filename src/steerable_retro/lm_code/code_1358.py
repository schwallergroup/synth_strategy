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
    This function detects a synthetic strategy involving amide coupling with
    halogenated aromatic compounds (specifically difluorobenzene).
    """
    # Track if we've found the pattern
    found_halogenated_aromatic_amide = False

    def dfs_traverse(node, depth=0):
        nonlocal found_halogenated_aromatic_amide

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is an amide coupling reaction - check multiple reaction types
                is_amide_coupling = (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction("Acylation of primary amines", rsmi)
                    or checker.check_reaction("Acylation of secondary amines", rsmi)
                    or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                )

                if is_amide_coupling:
                    print(f"Found amide coupling reaction at depth {depth}")

                    # Check for halogenated aromatic carboxylic acid in reactants
                    for reactant in reactants:
                        if checker.check_fg("Carboxylic acid", reactant) and "F" in reactant:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                # Check if it's a difluorobenzene
                                smiles = Chem.MolToSmiles(mol)

                                # Look for patterns indicating difluorobenzene
                                has_difluoro = False
                                for pattern in [
                                    "c(F)c(F)",
                                    "c(F)cc(F)",
                                    "c(F)ccc(F)",
                                    "c1c(F)cccc1F",
                                ]:
                                    if pattern in smiles:
                                        has_difluoro = True
                                        break

                                # Alternative check using substructure matching
                                if not has_difluoro:
                                    difluoro_patterns = [
                                        "c1c(F)c(F)ccc1",
                                        "c1c(F)cc(F)cc1",
                                        "c1c(F)ccc(F)c1",
                                    ]
                                    for pattern in difluoro_patterns:
                                        pattern_mol = Chem.MolFromSmiles(pattern)
                                        if pattern_mol and mol.HasSubstructMatch(pattern_mol):
                                            has_difluoro = True
                                            break

                                if has_difluoro:
                                    print(
                                        f"Found difluorinated aromatic carboxylic acid at depth {depth}"
                                    )

                                    # Check if product has an amide group and preserves the fluorine
                                    if (
                                        checker.check_fg("Primary amide", product)
                                        or checker.check_fg("Secondary amide", product)
                                        or checker.check_fg("Tertiary amide", product)
                                    ) and "F" in product:
                                        print(
                                            f"Confirmed halogenated aromatic amide coupling at depth {depth}"
                                        )
                                        found_halogenated_aromatic_amide = True
                                        return  # Found what we're looking for, no need to continue
                else:
                    # Check if this is an ester hydrolysis to carboxylic acid
                    is_ester_hydrolysis = checker.check_reaction(
                        "Ester saponification (methyl deprotection)", rsmi
                    ) or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)

                    if (
                        is_ester_hydrolysis
                        and "F" in product
                        and checker.check_fg("Carboxylic acid", product)
                    ):
                        print(
                            f"Found ester hydrolysis to fluorinated carboxylic acid at depth {depth}"
                        )
                        # This is part of the strategy - preparing the fluorinated acid for coupling

                    # Check for amide formation with difluorobenzene carboxylic acid
                    if "F" in product and (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    ):
                        for reactant in reactants:
                            if checker.check_fg("Carboxylic acid", reactant) and "F" in reactant:
                                mol = Chem.MolFromSmiles(reactant)
                                if mol:
                                    # Check for difluorobenzene pattern
                                    difluoro_patterns = [
                                        "c1c(F)c(F)ccc1",
                                        "c1c(F)cc(F)cc1",
                                        "c1c(F)ccc(F)c1",
                                    ]
                                    for pattern in difluoro_patterns:
                                        pattern_mol = Chem.MolFromSmiles(pattern)
                                        if pattern_mol and mol.HasSubstructMatch(pattern_mol):
                                            print(
                                                f"Found difluorinated aromatic carboxylic acid coupling at depth {depth}"
                                            )
                                            found_halogenated_aromatic_amide = True
                                            return
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if the pattern is found
    return found_halogenated_aromatic_amide

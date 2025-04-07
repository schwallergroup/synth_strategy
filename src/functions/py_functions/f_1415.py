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
    This function detects if the synthetic route involves azide formation
    and subsequent transformations using the azide group.
    """
    azide_molecules = set()  # Track molecules containing azide groups
    azide_reactions = []  # Track reactions involving azides
    azide_formation_reactions = []  # Track reactions where azide is formed

    def dfs_traverse(node, depth=0):
        nonlocal azide_reactions, azide_formation_reactions, azide_molecules

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check if this molecule contains an azide group
            if checker.check_fg("Azide", mol_smiles):
                azide_molecules.add(mol_smiles)
                print(f"Molecule with azide detected: {mol_smiles}")

        elif node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an azide formation reaction
            reactants_have_azide = any(checker.check_fg("Azide", r) for r in reactants)
            product_has_azide = checker.check_fg("Azide", product)

            # Azide formation: no azide in reactants, but azide in product
            if not reactants_have_azide and product_has_azide:
                azide_formation_reactions.append(depth)
                print(f"Azide formation detected at depth {depth}")

                # Check for specific azide formation reactions
                if (
                    checker.check_reaction("Formation of Azides from halogens", rsmi)
                    or checker.check_reaction(
                        "Formation of Azides from boronic acids", rsmi
                    )
                    or checker.check_reaction("Alcohol to azide", rsmi)
                    or checker.check_reaction("Amine to azide", rsmi)
                ):
                    print(f"Specific azide formation reaction detected: {rsmi}")

            # Azide transformation: azide in reactants
            elif reactants_have_azide:
                azide_reactions.append(depth)
                print(f"Azide transformation at depth {depth}")

                # Check for specific azide transformation reactions
                if (
                    checker.check_reaction(
                        "Huisgen alkyne-azide 1,3 dipolar cycloaddition", rsmi
                    )
                    or checker.check_reaction("Huisgen 1,3 dipolar cycloaddition", rsmi)
                    or checker.check_reaction(
                        "Huisgen alkene-azide 1,3 dipolar cycloaddition", rsmi
                    )
                    or checker.check_reaction(
                        "Azide-nitrile click cycloaddition to tetrazole", rsmi
                    )
                    or checker.check_reaction(
                        "Azide-nitrile click cycloaddition to triazole", rsmi
                    )
                    or checker.check_reaction(
                        "Azide to amine reduction (Staudinger)", rsmi
                    )
                ):
                    print(f"Specific azide transformation reaction detected: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have evidence of an azide-mediated strategy
    # Case 1: Azide is formed and then used in subsequent reactions
    if azide_formation_reactions and azide_reactions:
        print(f"Azide-mediated strategy detected: formed and then used in reactions")
        return True

    # Case 2: Starting material contains azide and it's used in multiple reactions
    if len(azide_molecules) > 0 and len(azide_reactions) >= 2:
        print(
            f"Azide-mediated strategy detected: starting with azide and used in multiple reactions"
        )
        return True

    # Case 3: Multiple azide transformations are present (even without explicit formation)
    if len(azide_reactions) >= 2:
        print(f"Azide-mediated strategy detected: multiple azide transformations")
        return True

    return False

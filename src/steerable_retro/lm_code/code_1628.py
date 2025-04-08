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
    Detects if the synthesis involves nitration of an aromatic ring.
    """
    # Track if we found aromatic nitration
    found_aromatic_nitration = False

    def dfs_traverse(node):
        nonlocal found_aromatic_nitration

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # First check if this is a known aromatic nitration reaction
            if (
                checker.check_reaction("Aromatic nitration with HNO3", rsmi)
                or checker.check_reaction("Aromatic nitration with NO3 salt", rsmi)
                or checker.check_reaction("Aromatic nitration with NO2 salt", rsmi)
                or checker.check_reaction("Aromatic nitration with alkyl NO2", rsmi)
            ):
                print(f"Found aromatic nitration reaction: {rsmi}")
                found_aromatic_nitration = True

            # If not identified by reaction type, check by analyzing reactants and products
            else:
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Check if product contains nitro group
                if checker.check_fg("Nitro group", product_part):
                    product_mol = Chem.MolFromSmiles(product_part)

                    # Check if the nitro group is on an aromatic ring
                    if product_mol:
                        nitro_aromatic_pattern = Chem.MolFromSmarts("c-[N+](=[O])[O-]")
                        if product_mol.HasSubstructMatch(nitro_aromatic_pattern):
                            # Check reactants for nitrating agents and aromatic rings
                            reactants = reactants_part.split(".")
                            has_nitrating_agent = False
                            has_aromatic_without_nitro = False

                            for reactant in reactants:
                                # Check for nitrating agents
                                if (
                                    "HNO3" in reactant
                                    or "NO3" in reactant
                                    or "NO2" in reactant
                                    or checker.check_fg("Nitro group", reactant)
                                ):
                                    has_nitrating_agent = True

                                # Check for aromatic rings without nitro groups
                                r_mol = Chem.MolFromSmiles(reactant)
                                if r_mol and r_mol.GetNumAtoms() > 0:
                                    has_aromatic = any(
                                        atom.GetIsAromatic() for atom in r_mol.GetAtoms()
                                    )
                                    has_nitro = checker.check_fg("Nitro group", reactant)

                                    if has_aromatic and not has_nitro:
                                        has_aromatic_without_nitro = True

                            if has_nitrating_agent and has_aromatic_without_nitro:
                                print(f"Found aromatic nitration by analysis: {rsmi}")
                                found_aromatic_nitration = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return found_aromatic_nitration

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
    This function detects if the synthetic route contains multiple cross-coupling reactions
    (Suzuki, Stille, etc.) for C-C bond formation between aromatic rings.
    """
    cross_coupling_count = 0
    cross_coupling_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic acids OTf",
        "Suzuki coupling with sulfonic esters",
        "Suzuki coupling with boronic esters OTf",
        "Suzuki coupling with boronic esters",
        "Stille reaction_vinyl",
        "Stille reaction_aryl",
        "Stille reaction_benzyl",
        "Stille reaction_allyl",
        "Stille reaction_vinyl OTf",
        "Stille reaction_aryl OTf",
        "Stille reaction_benzyl OTf",
        "Stille reaction_allyl OTf",
        "Stille reaction_other",
        "Stille reaction_other OTf",
        "Negishi coupling",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
        "Aryllithium cross-coupling",
        "decarboxylative_coupling",
    ]

    def is_cross_coupling(rsmi):
        try:
            # Check if the reaction is any of the known cross-coupling types
            for rxn_type in cross_coupling_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    print(f"Found {rxn_type} reaction: {rsmi}")
                    return True

            # If no specific reaction type matched, check for characteristic patterns
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Create RDKit molecules from SMILES
            reactant_mols = []
            for r in reactants:
                mol = Chem.MolFromSmiles(r)
                if mol:
                    reactant_mols.append(mol)

            # Check for presence of boron or tin reagents (Suzuki or Stille)
            boron_pattern = Chem.MolFromSmarts("[c,C][B]([O,OH])[O,OH]")
            tin_pattern = Chem.MolFromSmarts("[c,C][Sn]")
            halogen_pattern = Chem.MolFromSmarts("[c,C][Cl,Br,I]")

            has_boron = any(
                mol.HasSubstructMatch(boron_pattern) for mol in reactant_mols if mol
            )
            has_tin = any(
                mol.HasSubstructMatch(tin_pattern) for mol in reactant_mols if mol
            )
            has_halogen = any(
                mol.HasSubstructMatch(halogen_pattern) for mol in reactant_mols if mol
            )

            # If we have both a halogen compound and either boron or tin, likely a cross-coupling
            return has_halogen and (has_boron or has_tin)
        except Exception as e:
            print(f"Error in is_cross_coupling: {e}")
            return False

    def dfs_traverse(node):
        nonlocal cross_coupling_count

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                if is_cross_coupling(rsmi):
                    cross_coupling_count += 1
                    print(f"Cross-coupling reaction count: {cross_coupling_count}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total cross-coupling reactions found: {cross_coupling_count}")
    return cross_coupling_count >= 2

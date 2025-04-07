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
    Detects if the synthetic route uses amine alkylation as a key fragment coupling strategy.
    """
    has_amine_alkylation = False

    def dfs_traverse(node):
        nonlocal has_amine_alkylation

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]

            # Check for amine alkylation reactions using the checker function
            amine_alkylation_reactions = [
                "N-alkylation of primary amines with alkyl halides",
                "N-alkylation of secondary amines with alkyl halides",
                "Alkylation of amines",
                "Methylation with MeI_primary",
                "Methylation with MeI_secondary",
                "Methylation with MeI_tertiary",
                "Eschweiler-Clarke Primary Amine Methylation",
                "Eschweiler-Clarke Secondary Amine Methylation",
                "Reductive methylation of primary amine with formaldehyde",
                "N-methylation",
            ]

            for reaction_type in amine_alkylation_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Found amine alkylation coupling: {reaction_type}")
                    has_amine_alkylation = True
                    break

            # If no specific reaction type matched, check for general pattern
            if not has_amine_alkylation:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactants contain amine and alkyl halide
                has_amine = False
                has_alkyl_halide = False

                for reactant in reactants:
                    if (
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                        or checker.check_fg("Tertiary amine", reactant)
                    ):
                        has_amine = True
                    if (
                        checker.check_fg("Primary halide", reactant)
                        or checker.check_fg("Secondary halide", reactant)
                        or checker.check_fg("Tertiary halide", reactant)
                    ):
                        has_alkyl_halide = True

                # If both amine and alkyl halide are present in reactants, check if product has a new C-N bond
                if has_amine and has_alkyl_halide:
                    try:
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            # Check if the product has more C-N bonds than the reactants combined
                            reactant_mols = [
                                Chem.MolFromSmiles(r)
                                for r in reactants
                                if Chem.MolFromSmiles(r)
                            ]

                            # Count C-N bonds in reactants
                            reactant_cn_bonds = sum(
                                len(mol.GetSubstructMatches(Chem.MolFromSmarts("C-N")))
                                for mol in reactant_mols
                                if mol
                            )

                            # Count C-N bonds in product
                            product_cn_bonds = len(
                                product_mol.GetSubstructMatches(
                                    Chem.MolFromSmarts("C-N")
                                )
                            )

                            if product_cn_bonds > reactant_cn_bonds:
                                print(
                                    "Found amine alkylation coupling (pattern-based detection)"
                                )
                                has_amine_alkylation = True
                    except Exception as e:
                        print(f"Error in pattern-based detection: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_amine_alkylation

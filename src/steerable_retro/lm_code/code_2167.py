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
    Detects if the synthesis route involves nucleophilic aromatic substitution (SNAr)
    with a halogenated heterocycle.
    """
    found_snar = False

    def dfs_traverse(node, depth=0):
        nonlocal found_snar

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a known nucleophilic aromatic substitution reaction
            is_snar_reaction = (
                checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
            )

            # If it's a known SNAr reaction, we're done
            if is_snar_reaction:
                found_snar = True
                print(f"Found nucleophilic aromatic substitution reaction at depth: {depth}")
                return

            # If not a known reaction type, check for characteristic patterns
            has_halogenated_heterocycle = False
            has_nucleophile = False
            heterocycle_rings = [
                "pyridine",
                "pyrimidine",
                "pyrazine",
                "pyridazine",
                "furan",
                "thiophene",
                "oxazole",
                "thiazole",
                "isoxazole",
                "isothiazole",
                "imidazole",
                "pyrazole",
                "triazole",
                "tetrazole",
            ]

            for reactant in reactants:
                # Check for aromatic halide
                if checker.check_fg("Aromatic halide", reactant):
                    # Check if it's part of a heterocycle
                    for ring in heterocycle_rings:
                        if checker.check_ring(ring, reactant):
                            has_halogenated_heterocycle = True
                            print(f"Found halogenated heterocycle ({ring}) in reactant: {reactant}")
                            break

                # Check for nucleophiles
                if (
                    checker.check_fg("Phenol", reactant)
                    or checker.check_fg("Primary amine", reactant)
                    or checker.check_fg("Secondary amine", reactant)
                    or checker.check_fg("Aromatic thiol", reactant)
                    or checker.check_fg("Aliphatic thiol", reactant)
                ):
                    has_nucleophile = True
                    print(f"Found nucleophile in reactant: {reactant}")

            # Check if product has new C-O, C-N, or C-S bond
            product_has_c_hetero_bond = (
                checker.check_fg("Ether", product)
                or checker.check_fg("Aniline", product)
                or checker.check_fg("Secondary amine", product)
                or checker.check_fg("Tertiary amine", product)
                or checker.check_fg("Monosulfide", product)
            )

            # Check for electron-withdrawing groups that facilitate SNAr
            has_ewg = False
            for reactant in reactants:
                if (
                    checker.check_fg("Nitro group", reactant)
                    or checker.check_fg("Nitrile", reactant)
                    or checker.check_fg("Ketone", reactant)
                    or checker.check_fg("Ester", reactant)
                    or checker.check_fg("Carboxylic acid", reactant)
                ):
                    has_ewg = True
                    print(f"Found electron-withdrawing group in reactant: {reactant}")
                    break

            # Determine if this is an SNAr reaction based on all criteria
            if has_halogenated_heterocycle and has_nucleophile and product_has_c_hetero_bond:
                # Heterocycles are inherently electron-deficient, so EWG is not strictly required
                found_snar = True
                print(f"Found nucleophilic aromatic substitution pattern at depth: {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return found_snar

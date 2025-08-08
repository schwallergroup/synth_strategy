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
    Detects a strategy involving late-stage coupling of a heterocycle (benzimidazole)
    after earlier protection/deprotection and functional group manipulations.
    """
    # Track key features
    has_late_benzimidazole_coupling = False
    has_boc_protection = False
    has_nitrile_intermediate = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_benzimidazole_coupling, has_boc_protection, has_nitrile_intermediate

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for benzimidazole coupling or formation at late stage (depth 0-1)
            if depth <= 1:
                # Check if any reactant contains benzimidazole
                for reactant in reactants_smiles:
                    if checker.check_ring("benzimidazole", reactant):
                        has_late_benzimidazole_coupling = True
                        print(f"Found benzimidazole ring in reactant at depth {depth}")

                # Check if product contains benzimidazole (could be formed in the reaction)
                if checker.check_ring("benzimidazole", product_smiles):
                    has_late_benzimidazole_coupling = True
                    print(f"Found benzimidazole ring in product at depth {depth}")

                # Check for benzimidazole formation reactions
                if (
                    checker.check_reaction("benzimidazole_formation from aldehyde", rsmi)
                    or checker.check_reaction("benzimidazole_formation from acyl halide", rsmi)
                    or checker.check_reaction(
                        "benzimidazole_formation from ester/carboxylic acid", rsmi
                    )
                    or checker.check_reaction(
                        "{benzimidazole_derivatives_carboxylic-acid/ester}", rsmi
                    )
                    or checker.check_reaction("{benzimidazole_derivatives_aldehyde}", rsmi)
                ):
                    has_late_benzimidazole_coupling = True
                    print(f"Found benzimidazole formation reaction at depth {depth}")

            # Check for Boc protection/deprotection in early stages (depth > 1)
            if depth > 1:
                # Check for Boc protection reactions
                if (
                    checker.check_reaction("Boc amine protection", rsmi)
                    or checker.check_reaction("Boc amine protection explicit", rsmi)
                    or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                    or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                    or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                    or checker.check_reaction("Boc amine protection of primary amine", rsmi)
                ):
                    has_boc_protection = True
                    print(f"Found Boc protection reaction at depth {depth}")

                # Check for Boc deprotection reactions
                if (
                    checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                    or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                ):
                    has_boc_protection = True
                    print(f"Found Boc deprotection reaction at depth {depth}")

                # Check for Boc functional group in molecules
                for smiles in reactants_smiles + [product_smiles]:
                    if checker.check_fg("Boc", smiles):
                        has_boc_protection = True
                        print(f"Found Boc group in molecule at depth {depth}")

            # Check for nitrile intermediate
            for smiles in reactants_smiles + [product_smiles]:
                if checker.check_fg("Nitrile", smiles):
                    has_nitrile_intermediate = True
                    print(f"Found nitrile at depth {depth}")

        # For molecule nodes, check for benzimidazole and Boc groups
        elif node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for benzimidazole in late-stage molecules
            if depth <= 1 and checker.check_ring("benzimidazole", mol_smiles):
                has_late_benzimidazole_coupling = True
                print(f"Found benzimidazole ring in molecule at depth {depth}")

            # Check for Boc in early-stage molecules
            if depth > 1 and checker.check_fg("Boc", mol_smiles):
                has_boc_protection = True
                print(f"Found Boc group in molecule node at depth {depth}")

            # Check for nitrile
            if checker.check_fg("Nitrile", mol_smiles):
                has_nitrile_intermediate = True
                print(f"Found nitrile in molecule at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have late-stage heterocycle coupling and earlier protection
    strategy_present = has_late_benzimidazole_coupling and has_boc_protection

    # Additional confidence if we also have nitrile intermediate
    if strategy_present and has_nitrile_intermediate:
        print(
            "Confirmed full late-stage heterocycle coupling strategy with protection and nitrile intermediate"
        )
    elif strategy_present:
        print("Detected late-stage heterocycle coupling strategy with protection")
    else:
        if has_late_benzimidazole_coupling:
            print("Found late-stage benzimidazole coupling but no Boc protection")
        if has_boc_protection:
            print("Found Boc protection but no late-stage benzimidazole coupling")

    return strategy_present

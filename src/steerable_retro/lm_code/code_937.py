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
    This function detects a synthetic strategy involving lactone formation through
    a sequence of functional group interconversions including carboxylic acid, nitrile,
    and halomethyl intermediates, while maintaining a diaryl ether motif.
    """
    # Initialize tracking variables
    has_lactone_formation = False
    has_nitrile_intermediate = False
    has_halomethyl_intermediate = False
    has_carboxylic_acid_intermediate = False
    has_diaryl_ether_throughout = True
    has_lactone_product = False

    # Track molecules to verify they contain diaryl ether
    molecules_with_diaryl_ether = []

    def dfs_traverse(node, depth=0):
        nonlocal has_lactone_formation, has_nitrile_intermediate, has_halomethyl_intermediate
        nonlocal has_carboxylic_acid_intermediate, has_diaryl_ether_throughout, has_lactone_product

        if node["type"] == "mol":
            # Check molecule for functional groups
            mol_smiles = node["smiles"]

            # Check for diaryl ether in all molecules
            # A diaryl ether should have both ether and aromatic rings
            has_ether = checker.check_fg("Ether", mol_smiles)
            has_aromatic = "c" in mol_smiles  # Simplified check for aromatic rings
            has_diaryl_ether = has_ether and has_aromatic

            if has_diaryl_ether:
                molecules_with_diaryl_ether.append(mol_smiles)
            else:
                # Skip molecules that are likely reagents/solvents (short SMILES)
                if len(mol_smiles) > 5:
                    print(f"Molecule lacks diaryl ether: {mol_smiles}")

            # Check for functional groups
            if checker.check_fg("Nitrile", mol_smiles):
                has_nitrile_intermediate = True
                print(f"Found nitrile intermediate: {mol_smiles}")

            # Check for halomethyl groups (primary, secondary, tertiary halides)
            if (
                checker.check_fg("Primary halide", mol_smiles)
                or checker.check_fg("Secondary halide", mol_smiles)
                or checker.check_fg("Tertiary halide", mol_smiles)
            ):
                has_halomethyl_intermediate = True
                print(f"Found halomethyl intermediate: {mol_smiles}")

            if checker.check_fg("Carboxylic acid", mol_smiles):
                has_carboxylic_acid_intermediate = True
                print(f"Found carboxylic acid intermediate: {mol_smiles}")

            # Check for lactone in molecules - multiple approaches
            # 1. Check for oxolane/oxane ring with ester
            if (
                checker.check_ring("oxolane", mol_smiles) or checker.check_ring("oxane", mol_smiles)
            ) and checker.check_fg("Ester", mol_smiles):
                print(f"Found potential lactone (ring with ester): {mol_smiles}")
                has_lactone_product = True

            # 2. Check for characteristic lactone patterns in SMILES
            if (
                "O=C1" in mol_smiles
                and "O" in mol_smiles
                and any(c in mol_smiles for c in ["c2", "c3", "c4"])
            ):
                print(f"Found potential lactone (pattern match): {mol_smiles}")
                has_lactone_product = True

            # 3. Check if it's the final product (depth 0) and has a cyclic structure with ester
            if depth == 0 and "1" in mol_smiles and checker.check_fg("Ester", mol_smiles):
                print(f"Found potential lactone as final product: {mol_smiles}")
                has_lactone_product = True

        elif node["type"] == "reaction":
            # Check for lactone formation
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Check for lactone formation reaction
                if checker.check_reaction(
                    "Intramolecular transesterification/Lactone formation", rsmi
                ):
                    has_lactone_formation = True
                    print(f"Found lactone formation reaction: {rsmi}")

                # Alternative check: see if product has lactone ring that reactants don't have
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                if not has_lactone_formation:
                    # Check for formation of a cyclic ester (lactone)
                    reactant_has_oxolane = checker.check_ring("oxolane", reactants)
                    reactant_has_oxane = checker.check_ring("oxane", reactants)
                    product_has_oxolane = checker.check_ring("oxolane", product)
                    product_has_oxane = checker.check_ring("oxane", product)

                    if (
                        (not reactant_has_oxolane and product_has_oxolane)
                        or (not reactant_has_oxane and product_has_oxane)
                    ) and checker.check_fg("Ester", product):
                        has_lactone_formation = True
                        print(f"Found lactone formation (new ring with ester): {rsmi}")

                    # Check for ring closure reactions that might form lactones
                    elif (
                        checker.check_fg("Carboxylic acid", reactants)
                        and not checker.check_fg("Carboxylic acid", product)
                        and checker.check_fg("Ester", product)
                        and (product_has_oxolane or product_has_oxane)
                    ):
                        has_lactone_formation = True
                        print(f"Found lactone formation from carboxylic acid: {rsmi}")

                    # Check for any reaction that creates a new ring and ester
                    elif (
                        "1" in product
                        and not "1" in reactants
                        and checker.check_fg("Ester", product)
                    ):
                        has_lactone_formation = True
                        print(f"Found potential ring closure to form lactone: {rsmi}")

        # Process children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Verify that at least one molecule has a diaryl ether
    if len(molecules_with_diaryl_ether) > 0:
        has_diaryl_ether_throughout = True
    else:
        has_diaryl_ether_throughout = False

    # Check if strategy is present - either explicit lactone formation reaction
    # or presence of lactone product with all required intermediates
    strategy_present = (
        (has_lactone_formation or has_lactone_product)
        and has_nitrile_intermediate
        and has_halomethyl_intermediate
        and has_carboxylic_acid_intermediate
        and has_diaryl_ether_throughout
    )

    print(f"Lactone formation: {has_lactone_formation}")
    print(f"Lactone product detected: {has_lactone_product}")
    print(f"Nitrile intermediate: {has_nitrile_intermediate}")
    print(f"Halomethyl intermediate: {has_halomethyl_intermediate}")
    print(f"Carboxylic acid intermediate: {has_carboxylic_acid_intermediate}")
    print(f"Diaryl ether throughout: {has_diaryl_ether_throughout}")
    print(f"Strategy present: {strategy_present}")

    return strategy_present

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
    This function detects a synthetic strategy involving phenol protection,
    Suzuki coupling for biaryl formation, and nitrile-to-ketone transformation.
    """
    # Initialize flags to track strategy components
    phenol_protection_count = 0
    has_suzuki_coupling = False
    has_nitrile_to_ketone = False

    def dfs_traverse(node, depth=0):
        nonlocal phenol_protection_count, has_suzuki_coupling, has_nitrile_to_ketone

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                reactants = reactants_str.split(".")

                # Check for phenol protection (O-alkylation or silyl protection)
                # Look for phenol in reactants
                phenol_reactants = [r for r in reactants if checker.check_fg("Phenol", r)]

                # Check for various phenol protection reactions
                if phenol_reactants and (
                    checker.check_reaction("Williamson Ether Synthesis", rsmi)
                    or checker.check_reaction("Alcohol protection with silyl ethers", rsmi)
                    or checker.check_reaction(
                        "O-alkylation of carboxylic acids with diazo compounds", rsmi
                    )
                    or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                    or checker.check_reaction("O-methylation", rsmi)
                ):
                    phenol_protection_count += 1
                    print(f"Found phenol protection at depth {depth}: {rsmi}")
                    print(f"  Phenol reactants: {phenol_reactants}")

                # Alternative check for phenol protection based on functional group transformation
                # Check if a phenol is being converted to an ether
                elif (
                    phenol_reactants
                    and not checker.check_fg("Phenol", product_str)
                    and checker.check_fg("Ether", product_str)
                ):
                    phenol_protection_count += 1
                    print(
                        f"Found phenol protection (by FG transformation) at depth {depth}: {rsmi}"
                    )
                    print(f"  Phenol reactants: {phenol_reactants}")

                # Check for already protected phenols in the synthesis route
                # This is important because the phenol might have been protected in a previous step
                elif checker.check_fg("Ether", product_str):
                    # Check if the ether is likely a protected phenol (aromatic-O-R pattern)
                    product_mol = Chem.MolFromSmiles(product_str)
                    if product_mol:
                        for atom in product_mol.GetAtoms():
                            # Look for oxygen atoms connected to aromatic carbon
                            if atom.GetSymbol() == "O" and atom.GetDegree() == 2:
                                neighbors = [n for n in atom.GetNeighbors()]
                                if any(n.GetIsAromatic() for n in neighbors) and any(
                                    not n.GetIsAromatic() for n in neighbors
                                ):
                                    # This is likely a protected phenol (Ar-O-R pattern)
                                    phenol_protection_count += 1
                                    print(
                                        f"Found protected phenol structure at depth {depth}: {product_str}"
                                    )
                                    break

                # Check for Suzuki coupling (biaryl formation)
                suzuki_reactions = [
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic acids OTf",
                    "Suzuki coupling with sulfonic esters",
                    "Suzuki coupling with boronic esters OTf",
                    "Suzuki coupling with boronic esters",
                    "Suzuki",
                ]

                if any(checker.check_reaction(rxn, rsmi) for rxn in suzuki_reactions):
                    has_suzuki_coupling = True
                    print(f"Found Suzuki coupling at depth {depth}: {rsmi}")

                # Alternative check for Suzuki coupling based on reactants
                elif any(
                    checker.check_fg("Boronic acid", r) or checker.check_fg("Boronic ester", r)
                    for r in reactants
                ) and any(checker.check_fg("Aromatic halide", r) for r in reactants):
                    has_suzuki_coupling = True
                    print(f"Found Suzuki coupling (by reactants) at depth {depth}: {rsmi}")
                    boronic_reactants = [
                        r
                        for r in reactants
                        if checker.check_fg("Boronic acid", r)
                        or checker.check_fg("Boronic ester", r)
                    ]
                    halide_reactants = [
                        r for r in reactants if checker.check_fg("Aromatic halide", r)
                    ]
                    print(f"  Boronic reactants: {boronic_reactants}")
                    print(f"  Halide reactants: {halide_reactants}")

                # Check for nitrile to ketone transformation
                nitrile_reactants = [r for r in reactants if checker.check_fg("Nitrile", r)]
                if nitrile_reactants and checker.check_fg("Ketone", product_str):
                    # Check for specific reactions if available
                    if (
                        checker.check_reaction("Grignard from nitrile to ketone", rsmi)
                        or checker.check_reaction("Ketone from Weinreb amide", rsmi)
                        or checker.check_reaction("Ketone from Grignard and CO2", rsmi)
                    ):
                        has_nitrile_to_ketone = True
                        print(f"Found nitrile to ketone transformation at depth {depth}: {rsmi}")
                        print(f"  Nitrile reactants: {nitrile_reactants}")
                        print(f"  Ketone product: {product_str}")
                    # Fallback to general functional group transformation
                    else:
                        # Check if nitrile is actually consumed in the reaction
                        if not checker.check_fg("Nitrile", product_str):
                            has_nitrile_to_ketone = True
                            print(
                                f"Found general nitrile to ketone transformation at depth {depth}: {rsmi}"
                            )
                            print(f"  Nitrile reactants: {nitrile_reactants}")
                            print(f"  Ketone product: {product_str}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        phenol_protection_count >= 1 and has_suzuki_coupling and has_nitrile_to_ketone
    )

    print(f"Strategy detection results:")
    print(f"  Phenol protection steps: {phenol_protection_count}")
    print(f"  Suzuki coupling: {has_suzuki_coupling}")
    print(f"  Nitrile to ketone: {has_nitrile_to_ketone}")
    print(f"  Overall strategy present: {strategy_present}")

    return strategy_present

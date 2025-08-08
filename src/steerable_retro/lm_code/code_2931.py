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
    This function detects if the synthesis involves an ether bond (C-O-C)
    formation/disconnection as a key step.
    """
    found_ether_disconnection = False

    def has_alcohol(smiles):
        """Helper function to check if a molecule has any type of alcohol"""
        return (
            checker.check_fg("Primary alcohol", smiles)
            or checker.check_fg("Secondary alcohol", smiles)
            or checker.check_fg("Tertiary alcohol", smiles)
            or checker.check_fg("Aromatic alcohol", smiles)
            or checker.check_fg("Phenol", smiles)
        )

    def dfs_traverse(node, depth=0):
        nonlocal found_ether_disconnection

        print(f"Traversing node at depth {depth}: {node.get('type', 'unknown')}")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction SMILES: {rsmi}")

            # Check for Williamson ether synthesis or similar reactions
            if (
                checker.check_reaction("Williamson Ether Synthesis", rsmi)
                or checker.check_reaction("Williamson Ether Synthesis (intra to epoxy)", rsmi)
                or checker.check_reaction("{Williamson ether}", rsmi)
            ):
                found_ether_disconnection = True
                print(f"Found Williamson ether synthesis: {rsmi}")

            # Check for other ether formation reactions
            if (
                checker.check_reaction("Mitsunobu aryl ether", rsmi)
                or checker.check_reaction("Mitsunobu aryl ether (intramolecular)", rsmi)
                or checker.check_reaction("Alcohol to ether", rsmi)
                or checker.check_reaction(
                    "O-alkylation of carboxylic acids with diazo compounds", rsmi
                )
                or checker.check_reaction("O-alkylation of amides with diazo compounds", rsmi)
                or checker.check_reaction("Chan-Lam etherification", rsmi)
            ):
                found_ether_disconnection = True
                print(f"Found ether formation reaction: {rsmi}")

            # Check for ether cleavage reactions
            if (
                checker.check_reaction("Ether cleavage to primary alcohol", rsmi)
                or checker.check_reaction("Cleavage of methoxy ethers to alcohols", rsmi)
                or checker.check_reaction("Cleavage of alkoxy ethers to alcohols", rsmi)
            ):
                found_ether_disconnection = True
                print(f"Found ether cleavage reaction: {rsmi}")

            # Check for TMS ether protection/deprotection
            if (
                checker.check_reaction("Alcohol protection with silyl ethers", rsmi)
                or checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                or checker.check_reaction("Alcohol deprotection from silyl ethers (double)", rsmi)
                or checker.check_reaction("Alcohol deprotection from silyl ethers (diol)", rsmi)
                or checker.check_reaction("TMS deprotection from alkyne", rsmi)
            ):
                found_ether_disconnection = True
                print(f"Found silyl ether protection/deprotection: {rsmi}")

            # If no specific reaction is found, check for structural changes
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Reactants: {reactants}")
                print(f"Product: {product}")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if all(reactant_mols) and product_mol:
                    # Check if product has ether
                    product_has_ether = checker.check_fg("Ether", product)

                    # Check if reactants have ether
                    reactants_have_ether = any(checker.check_fg("Ether", r) for r in reactants)

                    # Check if reactants have alcohol
                    reactants_have_alcohol = any(has_alcohol(r) for r in reactants)

                    # Check if product has alcohol
                    product_has_alcohol = has_alcohol(product)

                    # Check for halides in reactants (potential alkylating agents)
                    reactants_have_halide = any(
                        checker.check_fg("Primary halide", r)
                        or checker.check_fg("Secondary halide", r)
                        or checker.check_fg("Tertiary halide", r)
                        or checker.check_fg("Aromatic halide", r)
                        for r in reactants
                    )

                    print(f"Product has ether: {product_has_ether}")
                    print(f"Reactants have ether: {reactants_have_ether}")
                    print(f"Reactants have alcohol: {reactants_have_alcohol}")
                    print(f"Product has alcohol: {product_has_alcohol}")
                    print(f"Reactants have halide: {reactants_have_halide}")

                    # Ether formation: alcohol + halide → ether
                    if product_has_ether and reactants_have_alcohol and reactants_have_halide:
                        found_ether_disconnection = True
                        print(f"Found ether formation (alcohol + halide): {rsmi}")

                    # Ether formation: alcohol → ether (any case where alcohols are consumed)
                    if product_has_ether and reactants_have_alcohol and not product_has_alcohol:
                        found_ether_disconnection = True
                        print(f"Found ether formation (alcohol consumed): {rsmi}")

                    # Ether disconnection: ether → alcohol
                    if reactants_have_ether and product_has_alcohol and not reactants_have_alcohol:
                        found_ether_disconnection = True
                        print(f"Found ether disconnection (ether to alcohol): {rsmi}")

                    # Check for benzofuran formation (common heterocyclic ether)
                    if checker.check_ring("benzofuran", product) and not any(
                        checker.check_ring("benzofuran", r) for r in reactants
                    ):
                        found_ether_disconnection = True
                        print(f"Found benzofuran formation (heterocyclic ether): {rsmi}")

                    # Check for furan formation
                    if checker.check_ring("furan", product) and not any(
                        checker.check_ring("furan", r) for r in reactants
                    ):
                        found_ether_disconnection = True
                        print(f"Found furan formation (heterocyclic ether): {rsmi}")

                    # Check for TMS ether protective group
                    product_has_tms_ether = checker.check_fg("TMS ether protective group", product)
                    reactants_have_tms_ether = any(
                        checker.check_fg("TMS ether protective group", r) for r in reactants
                    )

                    if (product_has_tms_ether and not reactants_have_tms_ether) or (
                        reactants_have_tms_ether and not product_has_tms_ether
                    ):
                        found_ether_disconnection = True
                        print(f"Found TMS ether formation/cleavage: {rsmi}")

                    # Check for silyl protective group (broader category)
                    product_has_silyl = checker.check_fg("Silyl protective group", product)
                    reactants_have_silyl = any(
                        checker.check_fg("Silyl protective group", r) for r in reactants
                    )

                    if (
                        product_has_silyl and not reactants_have_silyl and reactants_have_alcohol
                    ) or (reactants_have_silyl and not product_has_silyl and product_has_alcohol):
                        found_ether_disconnection = True
                        print(f"Found silyl protective group formation/cleavage: {rsmi}")

                    # Check for acetal/ketal formation (contains ether bonds)
                    product_has_acetal = checker.check_fg("Acetal/Ketal", product)
                    reactants_have_acetal = any(
                        checker.check_fg("Acetal/Ketal", r) for r in reactants
                    )

                    if (product_has_acetal and not reactants_have_acetal) or (
                        reactants_have_acetal and not product_has_acetal
                    ):
                        found_ether_disconnection = True
                        print(f"Found acetal/ketal formation/hydrolysis: {rsmi}")

            except Exception as e:
                print(f"Error analyzing reaction structure: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting traversal to find ether bond disconnection...")
    dfs_traverse(route)

    print(f"Ether bond disconnection found: {found_ether_disconnection}")
    return found_ether_disconnection

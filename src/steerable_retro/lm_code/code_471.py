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
    Detects if the route involves protecting group operations (protection and deprotection).
    Looks for various protection and deprotection reactions including Boc, silyl, and ether groups.
    """
    protection_found = False
    deprotection_found = False

    def dfs_traverse(node):
        nonlocal protection_found, deprotection_found

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for protection reactions
                protection_reactions = [
                    "Boc amine protection",
                    "Boc amine protection explicit",
                    "Boc amine protection with Boc anhydride",
                    "Boc amine protection (ethyl Boc)",
                    "Boc amine protection of secondary amine",
                    "Boc amine protection of primary amine",
                    "Alcohol protection with silyl ethers",
                    "Protection of carboxylic acid",
                    "Aldehyde or ketone acetalization",
                    "Diol acetalization",
                ]

                for rxn_name in protection_reactions:
                    if checker.check_reaction(rxn_name, rsmi):
                        print(f"Protection reaction detected: {rxn_name}")
                        protection_found = True
                        break

                # Check for deprotection reactions
                deprotection_reactions = [
                    "Boc amine deprotection",
                    "Boc amine deprotection of guanidine",
                    "Boc amine deprotection to NH-NH2",
                    "Alcohol deprotection from silyl ethers",
                    "Alcohol deprotection from silyl ethers (double)",
                    "Alcohol deprotection from silyl ethers (diol)",
                    "Cleavage of methoxy ethers to alcohols",
                    "Cleavage of alkoxy ethers to alcohols",
                    "Ether cleavage to primary alcohol",
                    "COOH ethyl deprotection",
                    "Hydroxyl benzyl deprotection",
                    "Carboxyl benzyl deprotection",
                    "Tert-butyl deprotection of amine",
                    "N-glutarimide deprotection",
                    "Phthalimide deprotection",
                    "TMS deprotection from alkyne",
                    "Acetal hydrolysis to diol",
                    "Acetal hydrolysis to aldehyde",
                    "Ketal hydrolysis to ketone",
                    "Deprotection of carboxylic acid",
                    "Ester saponification (methyl deprotection)",
                    "Ester saponification (alkyl deprotection)",
                    "Hydrogenolysis of amides/imides/carbamates",
                    "Hydrolysis of amides/imides/carbamates",
                ]

                for rxn_name in deprotection_reactions:
                    if checker.check_reaction(rxn_name, rsmi):
                        print(f"Deprotection reaction detected: {rxn_name}")
                        deprotection_found = True
                        break

                # If no specific reaction was found, check for functional group changes
                if not protection_found:
                    # Check for Carbamic ester protection (Cbz is a type of carbamic ester)
                    has_carbamic_ester_in_product = checker.check_fg("Carbamic ester", product)

                    if has_carbamic_ester_in_product:
                        has_carbamic_ester_in_reactants = any(
                            checker.check_fg("Carbamic ester", r) for r in reactants
                        )

                        if not has_carbamic_ester_in_reactants:
                            print("Carbamic ester protection detected through FG analysis")
                            protection_found = True

                    # Check for acetal/ketal protection
                    has_acetal_in_product = checker.check_fg("Acetal/Ketal", product)
                    if has_acetal_in_product:
                        has_acetal_in_reactants = any(
                            checker.check_fg("Acetal/Ketal", r) for r in reactants
                        )
                        has_aldehyde_or_ketone_in_reactants = any(
                            checker.check_fg("Aldehyde", r) or checker.check_fg("Ketone", r)
                            for r in reactants
                        )

                        if not has_acetal_in_reactants and has_aldehyde_or_ketone_in_reactants:
                            print("Acetal/ketal protection detected through FG analysis")
                            protection_found = True

                    # Check for ether protection
                    has_ether_in_product = checker.check_fg("Ether", product)
                    if has_ether_in_product:
                        has_ether_in_reactants = any(
                            checker.check_fg("Ether", r) for r in reactants
                        )
                        has_alcohol_in_reactants = any(
                            checker.check_fg("Primary alcohol", r)
                            or checker.check_fg("Secondary alcohol", r)
                            or checker.check_fg("Tertiary alcohol", r)
                            or checker.check_fg("Aromatic alcohol", r)
                            for r in reactants
                        )

                        if not has_ether_in_reactants and has_alcohol_in_reactants:
                            print("Ether protection detected through FG analysis")
                            protection_found = True

                if not deprotection_found:
                    # Check for ether deprotection by functional group analysis
                    product_has_alcohol = (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                        or checker.check_fg("Aromatic alcohol", product)
                    )

                    if product_has_alcohol:
                        # Check for ethers in reactants
                        reactants_have_ether = any(
                            checker.check_fg("TMS ether protective group", r)
                            or checker.check_fg("Silyl protective group", r)
                            or checker.check_fg("Ether", r)
                            for r in reactants
                        )

                        if reactants_have_ether and not any(
                            checker.check_fg("Primary alcohol", r)
                            or checker.check_fg("Secondary alcohol", r)
                            or checker.check_fg("Tertiary alcohol", r)
                            or checker.check_fg("Aromatic alcohol", r)
                            for r in reactants
                        ):
                            print("Ether deprotection detected through FG analysis")
                            deprotection_found = True

                    # Check for acetal/ketal deprotection
                    product_has_carbonyl = checker.check_fg(
                        "Aldehyde", product
                    ) or checker.check_fg("Ketone", product)
                    if product_has_carbonyl:
                        reactants_have_acetal = any(
                            checker.check_fg("Acetal/Ketal", r) for r in reactants
                        )
                        if reactants_have_acetal and not any(
                            checker.check_fg("Aldehyde", r) or checker.check_fg("Ketone", r)
                            for r in reactants
                        ):
                            print("Acetal/ketal deprotection detected through FG analysis")
                            deprotection_found = True

                    # Check for carboxylic acid deprotection from esters
                    product_has_carboxylic_acid = checker.check_fg("Carboxylic acid", product)
                    if product_has_carboxylic_acid:
                        reactants_have_ester = any(checker.check_fg("Ester", r) for r in reactants)
                        if reactants_have_ester and not any(
                            checker.check_fg("Carboxylic acid", r) for r in reactants
                        ):
                            print(
                                "Ester deprotection to carboxylic acid detected through FG analysis"
                            )
                            deprotection_found = True

                    # Check for amine deprotection from carbamates or amides
                    product_has_amine = (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                    )
                    if product_has_amine:
                        reactants_have_protected_amine = any(
                            checker.check_fg("Carbamic ester", r)
                            or checker.check_fg("Primary amide", r)
                            or checker.check_fg("Secondary amide", r)
                            or checker.check_fg("Tertiary amide", r)
                            or checker.check_fg("Boc", r)
                            for r in reactants
                        )
                        if reactants_have_protected_amine and not any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            or checker.check_fg("Tertiary amine", r)
                            for r in reactants
                        ):
                            print("Amine deprotection detected through FG analysis")
                            deprotection_found = True

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Protection found: {protection_found}, Deprotection found: {deprotection_found}")
    # Return True if both protection and deprotection are found
    return protection_found and deprotection_found

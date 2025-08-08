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
    This function detects if the synthetic route follows a linear synthesis strategy
    with sequential functional group transformations including esterification,
    reductive amination, and oxidation.
    """
    esterification_detected = False
    reductive_amination_detected = False
    oxidation_detected = False
    reaction_count = 0
    reaction_sequence = []

    def dfs_traverse(node, depth=0):
        nonlocal esterification_detected, reductive_amination_detected, oxidation_detected, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1
            # Add to reaction sequence (will reverse later to get forward synthesis order)
            reaction_sequence.append(node)

            # Get reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction {reaction_count}: {rsmi}")

            # Check for esterification
            if checker.check_reaction("Esterification of Carboxylic Acids", rsmi):
                # Verify carboxylic acid is converted to ester
                if any(
                    checker.check_fg("Carboxylic acid", r) for r in reactants
                ) and checker.check_fg("Ester", product):
                    esterification_detected = True
                    print(f"  Esterification detected in reaction {reaction_count}")

            # Check for reductive amination
            reductive_amination_rxns = [
                "Reductive amination with aldehyde",
                "Reductive amination with ketone",
                "Reductive amination with alcohol",
                "reductive amination",
            ]

            for rxn_type in reductive_amination_rxns:
                if checker.check_reaction(rxn_type, rsmi):
                    # Verify aldehyde/ketone and amine are converted
                    carbonyl_present = any(
                        checker.check_fg("Aldehyde", r) for r in reactants
                    ) or any(checker.check_fg("Ketone", r) for r in reactants)
                    amine_present = any(
                        checker.check_fg("Primary amine", r) for r in reactants
                    ) or any(checker.check_fg("Secondary amine", r) for r in reactants)
                    amine_product = (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                    )

                    if (
                        (carbonyl_present or any(checker.check_fg("Alcohol", r) for r in reactants))
                        and amine_present
                        and amine_product
                    ):
                        reductive_amination_detected = True
                        print(f"  Reductive amination detected in reaction {reaction_count}")
                        break

            # Check for oxidation reactions
            oxidation_reactions = [
                "Oxidation of aldehydes to carboxylic acids",
                "Oxidation of ketone to carboxylic acid",
                "Oxidation of alcohol to carboxylic acid",
                "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                "Oxidation of alkene to carboxylic acid",
                "Oxidation of alkene to aldehyde",
            ]

            # Check for oxidation - first try reaction type check
            oxidation_found = False
            for rxn_type in oxidation_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    oxidation_found = True
                    print(f"  Oxidation reaction type detected: {rxn_type}")
                    break

            # If reaction type check failed, try functional group analysis
            if not oxidation_found:
                # Check for aldehyde to carboxylic acid oxidation
                if any(checker.check_fg("Aldehyde", r) for r in reactants) and checker.check_fg(
                    "Carboxylic acid", product
                ):
                    oxidation_found = True
                    print(f"  Oxidation detected: Aldehyde to carboxylic acid")
                # Check for alcohol to aldehyde/ketone oxidation
                elif any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    for r in reactants
                ) and (
                    checker.check_fg("Aldehyde", product) or checker.check_fg("Ketone", product)
                ):
                    oxidation_found = True
                    print(f"  Oxidation detected: Alcohol to aldehyde/ketone")
                # Check for alcohol to carboxylic acid oxidation
                elif any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    for r in reactants
                ) and checker.check_fg("Carboxylic acid", product):
                    oxidation_found = True
                    print(f"  Oxidation detected: Alcohol to carboxylic acid")

            if oxidation_found:
                oxidation_detected = True
                print(f"  Oxidation confirmed in reaction {reaction_count}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Reverse reaction sequence to get forward synthesis order
    reaction_sequence.reverse()

    # Check if the synthesis is linear
    is_linear = True
    if len(reaction_sequence) >= 2:
        for i in range(len(reaction_sequence) - 1):
            current_reaction = reaction_sequence[i]["metadata"]["rsmi"]
            next_reaction = reaction_sequence[i + 1]["metadata"]["rsmi"]

            current_product = current_reaction.split(">")[-1]
            next_reactants = next_reaction.split(">")[0].split(".")

            print(f"Checking linearity between reactions {i+1} and {i+2}")
            print(f"  Product of reaction {i+1}: {current_product}")
            print(f"  Reactants of reaction {i+2}: {'.'.join(next_reactants)}")

            # Check if product of reaction i is a reactant in reaction i+1
            product_found = False
            for reactant in next_reactants:
                # Use atom mapping to check if product is part of the next reactant
                product_mol = Chem.MolFromSmiles(current_product)
                reactant_mol = Chem.MolFromSmiles(reactant)

                if product_mol and reactant_mol:
                    # Try to find common atom mappings between product and reactant
                    product_atoms = set()
                    for atom in product_mol.GetAtoms():
                        if atom.HasProp("molAtomMapNumber"):
                            product_atoms.add(atom.GetProp("molAtomMapNumber"))

                    reactant_atoms = set()
                    for atom in reactant_mol.GetAtoms():
                        if atom.HasProp("molAtomMapNumber"):
                            reactant_atoms.add(atom.GetProp("molAtomMapNumber"))

                    common_atoms = product_atoms.intersection(reactant_atoms)
                    if len(common_atoms) > 0:
                        product_found = True
                        print(f"  Found {len(common_atoms)} common mapped atoms")
                        break

                    # Fallback to substructure match for unmapped molecules
                    if len(product_atoms) == 0 or len(reactant_atoms) == 0:
                        if reactant_mol.HasSubstructMatch(
                            product_mol
                        ) or product_mol.HasSubstructMatch(reactant_mol):
                            product_found = True
                            print("  Found substructure match (fallback method)")
                            break

            if not product_found:
                is_linear = False
                print(f"Non-linear transition detected between reactions {i+1} and {i+2}")
                break
    else:
        is_linear = False  # Need at least 2 reactions to be linear

    # Check if it's a linear synthesis with multiple functional group transformations
    transformations_count = (
        esterification_detected + reductive_amination_detected + oxidation_detected
    )

    # For this function, we need at least 1 transformation in a linear synthesis with at least 2 reactions
    linear_with_transformations = reaction_count >= 2 and transformations_count >= 1 and is_linear

    print(
        f"Linear synthesis with functional group transformations detected: {linear_with_transformations}"
    )
    print(f"  - Total reactions: {reaction_count}")
    print(f"  - Is linear: {is_linear}")
    print(f"  - Transformations detected: {transformations_count}")
    print(f"  - Esterification: {esterification_detected}")
    print(f"  - Reductive amination: {reductive_amination_detected}")
    print(f"  - Oxidation: {oxidation_detected}")

    return linear_with_transformations

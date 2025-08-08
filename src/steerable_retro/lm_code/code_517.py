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
    Detects a synthesis strategy involving multiple functional group interconversions,
    particularly focusing on alcohol oxidations and halide transformations.
    """
    fg_interconversion_count = 0
    alcohol_oxidation_count = 0
    halide_transformation_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal fg_interconversion_count, alcohol_oxidation_count, halide_transformation_count

        if node["type"] == "reaction":
            # Extract reactants and product from reaction
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # In retrosynthesis, the product is the target and reactants are the precursors
                # But in reaction SMILES, reactants -> product is the forward direction

                # Check for alcohol oxidation reactions
                alcohol_oxidation_reactions = [
                    "Oxidation of aldehydes to carboxylic acids",
                    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                    "Oxidation of alcohol to carboxylic acid",
                    "Oxidation of ketone to carboxylic acid",
                    "Oxidation of primary alcohols",
                    "Oxidation of secondary alcohols",
                    "Oxidative esterification of primary alcohols",
                ]

                # Check for alcohol functional groups
                reactant_has_alcohol = any(
                    any(
                        checker.check_fg(fg, r)
                        for fg in [
                            "Primary alcohol",
                            "Secondary alcohol",
                            "Tertiary alcohol",
                            "Aromatic alcohol",
                            "Enol",
                        ]
                    )
                    for r in reactants_smiles
                )

                product_has_oxidized_group = (
                    checker.check_fg("Aldehyde", product_smiles)
                    or checker.check_fg("Ketone", product_smiles)
                    or checker.check_fg("Carboxylic acid", product_smiles)
                    or checker.check_fg("Ester", product_smiles)
                )

                # Check if this is an alcohol oxidation
                if (reactant_has_alcohol and product_has_oxidized_group) or any(
                    checker.check_reaction(rxn, rsmi) for rxn in alcohol_oxidation_reactions
                ):
                    alcohol_oxidation_count += 1
                    fg_interconversion_count += 1
                    print(f"Alcohol oxidation detected at depth {depth}, reaction: {rsmi}")

                # Check for halide transformations
                halide_fgs = [
                    "Primary halide",
                    "Secondary halide",
                    "Tertiary halide",
                    "Aromatic halide",
                    "Alkenyl halide",
                ]

                # Count halides in reactants and products
                reactant_halides = sum(
                    sum(1 for fg in halide_fgs if checker.check_fg(fg, r)) for r in reactants_smiles
                )

                product_halides = sum(
                    1 for fg in halide_fgs if checker.check_fg(fg, product_smiles)
                )

                # Check for halide transformation reactions
                halide_reactions = [
                    "Aromatic fluorination",
                    "Aromatic chlorination",
                    "Aromatic bromination",
                    "Aromatic iodination",
                    "Chlorination",
                    "Fluorination",
                    "Iodination",
                    "Bromination",
                    "Finkelstein reaction",
                    "Alcohol to chloride_SOCl2",
                    "Alcohol to chloride_HCl",
                    "Primary amine to fluoride",
                    "Primary amine to chloride",
                    "Primary amine to bromide",
                    "Primary amine to iodide",
                    "Appel reaction",
                    "Halodeboronation of boronic acids",
                    "Halodeboronation of boronic esters",
                    "Aromatic substitution of bromine by chlorine",
                    "Aromatic dehalogenation",
                    "Dehalogenation",
                ]

                # Check if this is a halide transformation
                if (reactant_halides != product_halides) or any(
                    checker.check_reaction(rxn, rsmi) for rxn in halide_reactions
                ):
                    halide_transformation_count += 1
                    fg_interconversion_count += 1
                    print(f"Halide transformation detected at depth {depth}, reaction: {rsmi}")

                # Check for other functional group interconversions
                other_fg_reactions = [
                    "Reduction of aldehydes and ketones to alcohols",
                    "Reduction of carboxylic acid to primary alcohol",
                    "Reduction of ester to primary alcohol",
                    "Reduction of nitrile to amine",
                    "Reduction of nitro groups to amines",
                    "Reduction of primary amides to amines",
                    "Reduction of secondary amides to amines",
                    "Reduction of tertiary amides to amines",
                    "Esterification of Carboxylic Acids",
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    "Ester saponification (methyl deprotection)",
                    "Ester saponification (alkyl deprotection)",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acylation of secondary amines with anhydrides",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Reductive amination with alcohol",
                    "Azide to amine reduction (Staudinger)",
                    "Formation of Azides from halogens",
                    "Amine to azide",
                    "Alcohol to azide",
                ]

                if any(checker.check_reaction(rxn, rsmi) for rxn in other_fg_reactions):
                    fg_interconversion_count += 1
                    print(
                        f"Other functional group interconversion detected at depth {depth}, reaction: {rsmi}"
                    )

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the strategy criteria are met - at least 3 total interconversions
    # with at least 1 alcohol oxidation and 1 halide transformation
    strategy_detected = (
        fg_interconversion_count >= 3
        and alcohol_oxidation_count >= 1
        and halide_transformation_count >= 1
    )

    print(f"Total functional group interconversions: {fg_interconversion_count}")
    print(f"Alcohol oxidations: {alcohol_oxidation_count}")
    print(f"Halide transformations: {halide_transformation_count}")

    if strategy_detected:
        print("Multiple functional group interconversion strategy detected")
    else:
        print(
            "Strategy not detected - requires at least 3 total interconversions with at least 1 alcohol oxidation and 1 halide transformation"
        )

    return strategy_detected

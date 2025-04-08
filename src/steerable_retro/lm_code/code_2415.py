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
    This function detects if the final step involves ester hydrolysis to form a carboxylic acid.
    """
    final_step_is_hydrolysis = False

    # First, identify the final product (target molecule)
    final_product = None
    if route["type"] == "mol" and not route.get("in_stock", False):
        final_product = route

    if final_product is None:
        print("Could not identify final product")
        return False

    print(f"Final product identified: {final_product['smiles']}")

    # Check if the final product has a carboxylic acid group
    has_carboxylic_acid = checker.check_fg("Carboxylic acid", final_product["smiles"])
    print(f"Final product has carboxylic acid: {has_carboxylic_acid}")

    if not has_carboxylic_acid:
        print("Final product does not have a carboxylic acid group")
        return False

    # Check the immediate precursor reaction (should be the only child of the final product)
    if "children" in final_product and len(final_product["children"]) > 0:
        final_reaction = final_product["children"][0]

        if final_reaction["type"] != "reaction":
            print("Child of final product is not a reaction")
            return False

        print("Examining final reaction step")

        # Extract reactants and product
        try:
            rsmi = final_reaction["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]
            print(f"Reaction SMILES: {rsmi}")
            print(f"Reactants: {reactants_smiles}")
            print(f"Product: {product_smiles}")

            # Check if this is an ester hydrolysis reaction - check multiple reaction types
            is_hydrolysis_reaction = checker.check_reaction(
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
            )
            is_saponification_methyl = checker.check_reaction(
                "Ester saponification (methyl deprotection)", rsmi
            )
            is_saponification_alkyl = checker.check_reaction(
                "Ester saponification (alkyl deprotection)", rsmi
            )
            is_cooh_deprotection = checker.check_reaction("COOH ethyl deprotection", rsmi)

            print(f"Is hydrolysis reaction: {is_hydrolysis_reaction}")
            print(f"Is methyl saponification: {is_saponification_methyl}")
            print(f"Is alkyl saponification: {is_saponification_alkyl}")
            print(f"Is COOH deprotection: {is_cooh_deprotection}")

            if (
                is_hydrolysis_reaction
                or is_saponification_methyl
                or is_saponification_alkyl
                or is_cooh_deprotection
            ):
                # Verify ester in reactants and carboxylic acid in product
                reactants_have_ester = any(
                    checker.check_fg("Ester", smi) for smi in reactants_smiles
                )
                product_has_acid = checker.check_fg("Carboxylic acid", product_smiles)

                print(f"Reactants have ester: {reactants_have_ester}")
                print(f"Product has carboxylic acid: {product_has_acid}")

                if reactants_have_ester and product_has_acid:
                    print("Late-stage ester hydrolysis detected!")
                    final_step_is_hydrolysis = True
            else:
                # Manual check for ester hydrolysis pattern if reaction type checks fail
                reactants_have_ester = any(
                    checker.check_fg("Ester", smi) for smi in reactants_smiles
                )
                product_has_acid = checker.check_fg("Carboxylic acid", product_smiles)

                print(f"Manual check - Reactants have ester: {reactants_have_ester}")
                print(f"Manual check - Product has carboxylic acid: {product_has_acid}")

                # Check if reactants contain common hydrolysis reagents
                hydrolysis_reagents = ["CCO", "O", "[Na+]", "[OH-]", "NaOH", "KOH", "LiOH"]
                has_hydrolysis_reagents = any(
                    reagent in rsmi.split(">")[1] for reagent in hydrolysis_reagents
                )

                if reactants_have_ester and product_has_acid and has_hydrolysis_reagents:
                    print("Late-stage ester hydrolysis detected through manual pattern check!")
                    final_step_is_hydrolysis = True
        except Exception as e:
            print(f"Error processing reaction: {e}")
    else:
        print("Final product has no precursor reactions")

    print(f"Final result: {final_step_is_hydrolysis}")
    return final_step_is_hydrolysis

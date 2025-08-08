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
    This function detects a synthetic strategy involving late-stage amide to ester conversion.
    Note: The route is traversed retrosynthetically, so we're looking for ester-to-amide
    conversions in the retrosynthetic direction.
    """
    amide_to_ester_conversion = False
    max_depth_for_late_stage = 2  # Define what "late-stage" means in terms of depth

    def dfs_traverse(node, depth=0):
        nonlocal amide_to_ester_conversion

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Only consider late-stage reactions
            if depth <= max_depth_for_late_stage:
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check for specific reactions that might convert amides to esters (in forward direction)
                # In retrosynthesis, these would be ester-to-amide conversions
                if checker.check_reaction("Schotten-Baumann to ester", rsmi):
                    print(f"Detected Schotten-Baumann to ester at depth {depth}")
                    # Need to verify this is actually converting an ester to amide in retrosynthesis
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # In retrosynthesis, check for ester in reactants and amide in product
                    reactant_has_ester = False
                    for r_smiles in reactants_smiles:
                        if checker.check_fg("Ester", r_smiles):
                            print(f"Found ester in reactant: {r_smiles}")
                            reactant_has_ester = True
                            break

                    # Check for amide in product
                    product_has_amide = (
                        checker.check_fg("Primary amide", product_smiles)
                        or checker.check_fg("Secondary amide", product_smiles)
                        or checker.check_fg("Tertiary amide", product_smiles)
                    )
                    if product_has_amide:
                        print(f"Found amide in product: {product_smiles}")

                    if reactant_has_ester and product_has_amide:
                        print(f"Confirmed amide to ester conversion strategy at depth {depth}")
                        amide_to_ester_conversion = True
                        return

                # Check for other relevant reactions
                if (
                    checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    or checker.check_reaction("Transesterification", rsmi)
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                        rsmi,
                    )
                ):
                    # These reactions could potentially be used for amide to ester conversion
                    # In retrosynthesis, need to verify ester is being converted to amide

                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # In retrosynthesis, check for ester in reactants
                    reactant_has_ester = False
                    for r_smiles in reactants_smiles:
                        if checker.check_fg("Ester", r_smiles):
                            print(f"Found ester in reactant: {r_smiles}")
                            reactant_has_ester = True
                            break

                    # Check for amide in product
                    product_has_amide = (
                        checker.check_fg("Primary amide", product_smiles)
                        or checker.check_fg("Secondary amide", product_smiles)
                        or checker.check_fg("Tertiary amide", product_smiles)
                    )
                    if product_has_amide:
                        print(f"Found amide in product: {product_smiles}")

                    # Check for alcohol in reactants (common reagent for esterification)
                    reactant_has_alcohol = False
                    for r_smiles in reactants_smiles:
                        if (
                            checker.check_fg("Primary alcohol", r_smiles)
                            or checker.check_fg("Secondary alcohol", r_smiles)
                            or checker.check_fg("Tertiary alcohol", r_smiles)
                            or checker.check_fg("Aromatic alcohol", r_smiles)
                        ):
                            print(f"Found alcohol in reactant: {r_smiles}")
                            reactant_has_alcohol = True
                            break

                    if reactant_has_ester and product_has_amide:
                        # Additional check to ensure we're not just adding an amide group elsewhere
                        ester_count_reactants = sum(
                            1 for r in reactants_smiles if checker.check_fg("Ester", r)
                        )

                        ester_count_product = 1 if checker.check_fg("Ester", product_smiles) else 0

                        print(
                            f"Ester count in reactants: {ester_count_reactants}, in product: {ester_count_product}"
                        )

                        # In retrosynthesis, if ester count decreased and amide is present in product,
                        # it's likely a conversion from ester to amide (which represents amide to ester in forward synthesis)
                        if ester_count_reactants > ester_count_product:
                            print(f"Detected amide to ester conversion strategy at depth {depth}")
                            amide_to_ester_conversion = True
                            return

                # If no specific reaction matched, check for structural changes
                if not amide_to_ester_conversion:
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # In retrosynthesis, check for ester in reactants
                    reactant_has_ester = False
                    for r_smiles in reactants_smiles:
                        if checker.check_fg("Ester", r_smiles):
                            print(f"Found ester in reactant: {r_smiles}")
                            reactant_has_ester = True
                            break

                    # Check for amide in product
                    product_has_amide = (
                        checker.check_fg("Primary amide", product_smiles)
                        or checker.check_fg("Secondary amide", product_smiles)
                        or checker.check_fg("Tertiary amide", product_smiles)
                    )
                    if product_has_amide:
                        print(f"Found amide in product: {product_smiles}")

                    if reactant_has_ester and product_has_amide:
                        # Additional check to ensure we're not just adding an amide group elsewhere
                        ester_count_reactants = sum(
                            1 for r in reactants_smiles if checker.check_fg("Ester", r)
                        )

                        ester_count_product = 1 if checker.check_fg("Ester", product_smiles) else 0

                        print(
                            f"Ester count in reactants: {ester_count_reactants}, in product: {ester_count_product}"
                        )

                        # In retrosynthesis, if ester count decreased and amide is present in product,
                        # it's likely a conversion from ester to amide (which represents amide to ester in forward synthesis)
                        if ester_count_reactants > ester_count_product:
                            print(f"Detected amide to ester conversion strategy at depth {depth}")
                            amide_to_ester_conversion = True
                            return

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: amide_to_ester_conversion = {amide_to_ester_conversion}")
    return amide_to_ester_conversion

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
    This function detects if the synthesis involves formation of an aldehyde from an alcohol.
    """
    aldehyde_formation_found = False

    def dfs_traverse(node):
        nonlocal aldehyde_formation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for specific reaction types directly (comprehensive list)
            if checker.check_reaction(
                "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
            ):
                print(f"Alcohol to aldehyde oxidation reaction detected: {rsmi}")
                aldehyde_formation_found = True
            elif checker.check_reaction("Alkene oxidation to aldehyde", rsmi):
                print(f"Alkene to aldehyde oxidation reaction detected: {rsmi}")
                aldehyde_formation_found = True
            elif checker.check_reaction("Hydration of alkyne to aldehyde", rsmi):
                print(f"Alkyne to aldehyde hydration reaction detected: {rsmi}")
                aldehyde_formation_found = True
            elif checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi):
                # This is the reverse direction (retrosynthetic perspective)
                print(f"Carboxylic acid to aldehyde reduction detected: {rsmi}")
                aldehyde_formation_found = True
            elif checker.check_reaction(
                "Reduction of ester to primary alcohol", rsmi
            ) and checker.check_fg("Aldehyde", product):
                # Ester reduction can sometimes lead to aldehyde formation
                print(f"Ester reduction to aldehyde detected: {rsmi}")
                aldehyde_formation_found = True
            else:
                # Check for aldehyde in product
                aldehyde_in_product = False
                try:
                    if checker.check_fg("Aldehyde", product):
                        print(f"Aldehyde found in product: {product}")
                        aldehyde_in_product = True
                    elif checker.check_fg("Formaldehyde", product):
                        print(f"Formaldehyde found in product: {product}")
                        aldehyde_in_product = True
                except Exception as e:
                    print(f"Error checking for aldehyde in product: {e}")

                # If aldehyde is in product, check reactants for potential precursors
                if aldehyde_in_product:
                    # Check for alcohols in reactants (potential oxidation)
                    for reactant in reactants:
                        try:
                            if checker.check_fg("Primary alcohol", reactant) or checker.check_fg(
                                "Secondary alcohol", reactant
                            ):
                                print(f"Alcohol found in reactant: {reactant}")
                                # This is likely an alcohol to aldehyde transformation
                                print(
                                    f"Potential alcohol to aldehyde transformation detected: {rsmi}"
                                )
                                aldehyde_formation_found = True
                                break
                        except Exception as e:
                            print(f"Error checking for alcohol: {e}")

                    # Check for other potential aldehyde-forming precursors
                    if not aldehyde_formation_found:
                        for reactant in reactants:
                            try:
                                if checker.check_fg("Alkene", reactant):
                                    print(
                                        f"Potential alkene to aldehyde transformation detected: {rsmi}"
                                    )
                                    aldehyde_formation_found = True
                                    break
                                elif checker.check_fg("Alkyne", reactant):
                                    print(
                                        f"Potential alkyne to aldehyde transformation detected: {rsmi}"
                                    )
                                    aldehyde_formation_found = True
                                    break
                                elif checker.check_fg("Nitrile", reactant):
                                    print(
                                        f"Potential nitrile to aldehyde transformation detected: {rsmi}"
                                    )
                                    aldehyde_formation_found = True
                                    break
                            except Exception as e:
                                print(f"Error checking reactant: {e}")

                # Check for retrosynthetic perspective (aldehyde in reactant becomes something else)
                for reactant in reactants:
                    try:
                        if checker.check_fg("Aldehyde", reactant) and not checker.check_fg(
                            "Aldehyde", product
                        ):
                            print(f"Aldehyde in reactant transformed in product: {rsmi}")
                            aldehyde_formation_found = True
                            break
                    except Exception as e:
                        print(f"Error checking for aldehyde in reactant: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return aldehyde_formation_found

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
    Detects if the synthesis route employs a sequence of thiol functional group manipulations,
    including protection/deprotection and oxidation steps.
    """
    # Track the functional group states we observe
    observed_states = {
        "alcohol": False,  # Primary, Secondary, or Tertiary alcohol
        "mesylate": False,  # Mesylate
        "thioacetate": False,  # Thioacetate
        "thiol": False,  # Aliphatic or Aromatic thiol
        "methylthioether": False,  # Monosulfide
        "sulfone": False,  # Sulfone
    }

    # Track if we've seen relevant thiol manipulation reactions
    thiol_reactions_observed = False

    # Track transformations between functional groups
    transformations = []

    def dfs_traverse(node, depth=0):
        nonlocal thiol_reactions_observed

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for thiol-related reactions
            thiol_reaction = False

            # Track functional groups in reactants and product
            reactant_fgs = set()
            product_fgs = set()

            # Check reactants for functional groups
            for r in reactants:
                # Check for alcohols
                if (
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                ):
                    observed_states["alcohol"] = True
                    reactant_fgs.add("alcohol")
                    print(f"Alcohol detected in reactant at depth {depth}")

                # Check for mesylate
                if checker.check_fg("Mesylate", r):
                    observed_states["mesylate"] = True
                    reactant_fgs.add("mesylate")
                    print(f"Mesylate detected in reactant at depth {depth}")

                # Check for thioacetate
                if "SC(=O)C" in r:
                    observed_states["thioacetate"] = True
                    reactant_fgs.add("thioacetate")
                    print(f"Thioacetate detected in reactant at depth {depth}")

                # Check for thiols
                if checker.check_fg("Aliphatic thiol", r) or checker.check_fg("Aromatic thiol", r):
                    observed_states["thiol"] = True
                    reactant_fgs.add("thiol")
                    print(f"Thiol detected in reactant at depth {depth}")

                # Check for methylthioether/monosulfide
                if checker.check_fg("Monosulfide", r):
                    observed_states["methylthioether"] = True
                    reactant_fgs.add("methylthioether")
                    print(f"Methylthioether detected in reactant at depth {depth}")

                # Check for sulfone
                if checker.check_fg("Sulfone", r):
                    observed_states["sulfone"] = True
                    reactant_fgs.add("sulfone")
                    print(f"Sulfone detected in reactant at depth {depth}")

            # Check product for functional groups
            # Check for alcohols
            if (
                checker.check_fg("Primary alcohol", product)
                or checker.check_fg("Secondary alcohol", product)
                or checker.check_fg("Tertiary alcohol", product)
            ):
                observed_states["alcohol"] = True
                product_fgs.add("alcohol")
                print(f"Alcohol detected in product at depth {depth}")

            # Check for mesylate
            if checker.check_fg("Mesylate", product):
                observed_states["mesylate"] = True
                product_fgs.add("mesylate")
                print(f"Mesylate detected in product at depth {depth}")

            # Check for thioacetate
            if "SC(=O)C" in product:
                observed_states["thioacetate"] = True
                product_fgs.add("thioacetate")
                print(f"Thioacetate detected in product at depth {depth}")

            # Check for thiols
            if checker.check_fg("Aliphatic thiol", product) or checker.check_fg(
                "Aromatic thiol", product
            ):
                observed_states["thiol"] = True
                product_fgs.add("thiol")
                print(f"Thiol detected in product at depth {depth}")

            # Check for methylthioether/monosulfide
            if checker.check_fg("Monosulfide", product):
                observed_states["methylthioether"] = True
                product_fgs.add("methylthioether")
                print(f"Methylthioether detected in product at depth {depth}")

            # Check for sulfone
            if checker.check_fg("Sulfone", product):
                observed_states["sulfone"] = True
                product_fgs.add("sulfone")
                print(f"Sulfone detected in product at depth {depth}")

            # Check for specific thiol-related reactions
            if (
                checker.check_reaction("S-alkylation of thiols", rsmi)
                or checker.check_reaction("S-alkylation of thiols (ethyl)", rsmi)
                or checker.check_reaction("S-alkylation of thiols with alcohols", rsmi)
                or checker.check_reaction("S-alkylation of thiols with alcohols (ethyl)", rsmi)
            ):
                thiol_reaction = True
                print(f"S-alkylation reaction detected: {rsmi}")

            # Sulfanyl to sulfinyl/sulfone reactions
            elif (
                checker.check_reaction("Sulfanyl to sulfinyl", rsmi)
                or checker.check_reaction("Sulfanyl to sulfinyl_peroxide", rsmi)
                or checker.check_reaction("Sulfanyl to sulfinyl_H2O2", rsmi)
                or checker.check_reaction("Sulfanyl to sulfinyl_SO3-", rsmi)
                or checker.check_reaction("Sulfanyl to sulfinyl_sulfonyl", rsmi)
            ):
                thiol_reaction = True
                print(f"Sulfanyl oxidation reaction detected: {rsmi}")

            # Mesylate formation
            elif checker.check_reaction("Formation of Sulfonic Esters", rsmi):
                thiol_reaction = True
                print(f"Mesylate formation detected: {rsmi}")

            # Detect transformations by comparing functional groups in reactants vs product
            # Since we're traversing retrosynthetically, the product is the starting material
            # and reactants are the synthetic targets

            # Alcohol to mesylate transformation
            if "alcohol" in reactant_fgs and "mesylate" in product_fgs:
                transformations.append(("alcohol", "mesylate"))
                thiol_reaction = True
                print(f"Transformation: alcohol → mesylate at depth {depth}")

            # Mesylate to thioacetate transformation
            if "mesylate" in reactant_fgs and "thioacetate" in product_fgs:
                transformations.append(("mesylate", "thioacetate"))
                thiol_reaction = True
                print(f"Transformation: mesylate → thioacetate at depth {depth}")

            # Thioacetate to thiol transformation
            if "thioacetate" in reactant_fgs and "thiol" in product_fgs:
                transformations.append(("thioacetate", "thiol"))
                thiol_reaction = True
                print(f"Transformation: thioacetate → thiol at depth {depth}")

            # Thiol to methylthioether transformation
            if "thiol" in reactant_fgs and "methylthioether" in product_fgs:
                transformations.append(("thiol", "methylthioether"))
                thiol_reaction = True
                print(f"Transformation: thiol → methylthioether at depth {depth}")

            # Methylthioether to sulfone transformation
            if "methylthioether" in reactant_fgs and "sulfone" in product_fgs:
                transformations.append(("methylthioether", "sulfone"))
                thiol_reaction = True
                print(f"Transformation: methylthioether → sulfone at depth {depth}")

            # Also check the reverse direction (since we're traversing retrosynthetically)
            # Mesylate to alcohol transformation
            if "mesylate" in reactant_fgs and "alcohol" in product_fgs:
                transformations.append(("mesylate", "alcohol"))
                thiol_reaction = True
                print(f"Transformation: mesylate → alcohol at depth {depth}")

            # Thioacetate to mesylate transformation
            if "thioacetate" in reactant_fgs and "mesylate" in product_fgs:
                transformations.append(("thioacetate", "mesylate"))
                thiol_reaction = True
                print(f"Transformation: thioacetate → mesylate at depth {depth}")

            # Thiol to thioacetate transformation
            if "thiol" in reactant_fgs and "thioacetate" in product_fgs:
                transformations.append(("thiol", "thioacetate"))
                thiol_reaction = True
                print(f"Transformation: thiol → thioacetate at depth {depth}")

            # Methylthioether to thiol transformation
            if "methylthioether" in reactant_fgs and "thiol" in product_fgs:
                transformations.append(("methylthioether", "thiol"))
                thiol_reaction = True
                print(f"Transformation: methylthioether → thiol at depth {depth}")

            # Sulfone to methylthioether transformation
            if "sulfone" in reactant_fgs and "methylthioether" in product_fgs:
                transformations.append(("sulfone", "methylthioether"))
                thiol_reaction = True
                print(f"Transformation: sulfone → methylthioether at depth {depth}")

            if thiol_reaction:
                thiol_reactions_observed = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Print summary of observed states
    print("Observed functional group states:")
    for state, observed in observed_states.items():
        print(f"  {state}: {observed}")
    print(f"Thiol manipulation reactions observed: {thiol_reactions_observed}")
    print(f"Transformations detected: {transformations}")

    # Check if we've observed at least 4 of the 6 states AND relevant reactions
    states_count = sum(1 for state in observed_states.values() if state)
    print(f"Total states observed: {states_count}/6")

    # Check for at least 2 transformations in the thiol manipulation pathway
    transformation_count = len(transformations)
    print(f"Total transformations observed: {transformation_count}")

    # Return true if we have enough states and transformations
    return states_count >= 4 and transformation_count >= 2

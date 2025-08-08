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
    This function detects a bidirectional functional group manipulation strategy
    involving an ester → acid → ester sequence.
    """
    # Track transformations in sequence with molecule information
    transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Depth {depth}, Reaction: {rsmi}")

                # Check for ester hydrolysis (ester → acid)
                if (
                    checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                    or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                    or checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                    )
                ):

                    # Check if reactants contain ester and product contains carboxylic acid
                    reactant_has_ester = any(checker.check_fg("Ester", r) for r in reactants_smiles)
                    product_has_acid = checker.check_fg("Carboxylic acid", product_smiles)

                    print(
                        f"Ester hydrolysis check - Reactant has ester: {reactant_has_ester}, Product has acid: {product_has_acid}"
                    )

                    if reactant_has_ester and product_has_acid:
                        print(f"Detected ester to acid transformation at depth {depth}")
                        transformations.append(("ester_to_acid", product_smiles, depth))

                # Check for esterification (acid → ester)
                if (
                    checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    or checker.check_reaction("Protection of carboxylic acid", rsmi)
                    or checker.check_reaction(
                        "O-alkylation of carboxylic acids with diazo compounds", rsmi
                    )
                ):

                    # Check if reactants contain carboxylic acid and product contains ester
                    reactant_has_acid = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                    )
                    product_has_ester = checker.check_fg("Ester", product_smiles)

                    print(
                        f"Esterification check - Reactant has acid: {reactant_has_acid}, Product has ester: {product_has_ester}"
                    )

                    if reactant_has_acid and product_has_ester:
                        print(f"Detected acid to ester transformation at depth {depth}")
                        transformations.append(("acid_to_ester", product_smiles, depth))

                # Additional check for any reaction that converts ester to acid
                if not transformations or transformations[-1][0] != "ester_to_acid":
                    reactant_has_ester = any(checker.check_fg("Ester", r) for r in reactants_smiles)
                    product_has_acid = checker.check_fg("Carboxylic acid", product_smiles)

                    if reactant_has_ester and product_has_acid:
                        print(f"Detected generic ester to acid transformation at depth {depth}")
                        transformations.append(("ester_to_acid", product_smiles, depth))

                # Additional check for any reaction that converts acid to ester
                if not transformations or transformations[-1][0] != "acid_to_ester":
                    reactant_has_acid = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                    )
                    product_has_ester = checker.check_fg("Ester", product_smiles)

                    if reactant_has_acid and product_has_ester:
                        print(f"Detected generic acid to ester transformation at depth {depth}")
                        transformations.append(("acid_to_ester", product_smiles, depth))

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Sort transformations by depth to ensure correct sequence analysis
    transformations.sort(key=lambda x: x[2])

    print(f"All transformations: {transformations}")

    # Check for the sequence ester→acid→ester in the synthetic route
    for i in range(len(transformations) - 1):
        if (
            transformations[i][0] == "ester_to_acid"
            and transformations[i + 1][0] == "acid_to_ester"
        ):
            print(
                f"Found ester→acid→ester sequence at depths {transformations[i][2]} and {transformations[i+1][2]}"
            )
            return True

    # Also check for the reverse sequence (acid→ester→acid)
    for i in range(len(transformations) - 1):
        if (
            transformations[i][0] == "acid_to_ester"
            and transformations[i + 1][0] == "ester_to_acid"
        ):
            print(
                f"Found acid→ester→acid sequence at depths {transformations[i][2]} and {transformations[i+1][2]}"
            )
            return True

    return False

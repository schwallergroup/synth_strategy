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
    This function detects if the synthetic route involves formation of an alcohol
    in the final step (depth 0).
    """
    late_stage_alcohol = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_alcohol

        print(f"Traversing node type: {node['type']} at depth: {depth}")

        # Check for reaction nodes at depth 1 (final reaction step)
        if node["type"] == "reaction" and depth == 1:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing final reaction: {rsmi}")

                # Check if product has alcohol group
                has_alcohol_in_product = False
                alcohol_types = [
                    "Primary alcohol",
                    "Secondary alcohol",
                    "Tertiary alcohol",
                    "Aromatic alcohol",
                    "Phenol",
                    "Enol",
                ]

                for alcohol_type in alcohol_types:
                    if checker.check_fg(alcohol_type, product_smiles):
                        print(f"Found {alcohol_type} in product: {product_smiles}")
                        has_alcohol_in_product = True
                        break

                if has_alcohol_in_product:
                    # Check if any reactant already has the alcohol group
                    has_alcohol_in_reactants = False
                    for reactant in reactants_smiles:
                        for alcohol_type in alcohol_types:
                            if checker.check_fg(alcohol_type, reactant):
                                print(f"Found {alcohol_type} in reactant: {reactant}")
                                has_alcohol_in_reactants = True
                                break
                        if has_alcohol_in_reactants:
                            break

                    # Check if this is an alcohol-forming reaction
                    alcohol_forming_reactions = [
                        "Reduction of aldehydes and ketones to alcohols",
                        "Reduction of ester to primary alcohol",
                        "Reduction of ketone to secondary alcohol",
                        "Reduction of carboxylic acid to primary alcohol",
                        "Primary alkyl halide to alcohol",
                        "Secondary alkyl halide to alcohol",
                        "anti-Markovnikov alkene hydration to alcohol",
                        "Markovnikov alkene hydration to alcohol",
                        "Grignard from aldehyde to alcohol",
                        "Grignard from ketone to alcohol",
                        "Alkene to diol",
                        "Oxidation of boronic acids",
                        "Oxidation of boronic esters",
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                    ]

                    is_alcohol_forming_reaction = False
                    for reaction_type in alcohol_forming_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(f"Detected alcohol-forming reaction: {reaction_type}")
                            is_alcohol_forming_reaction = True
                            break

                    # If not in our predefined list, check if alcohol is formed by comparing reactants and products
                    if not is_alcohol_forming_reaction and not has_alcohol_in_reactants:
                        print(
                            "Alcohol appears to be formed in this reaction (not in predefined list)"
                        )
                        is_alcohol_forming_reaction = True

                    # Alcohol formation confirmed if:
                    # 1. Product has alcohol
                    # 2. Either reactants don't have alcohol OR it's a known alcohol-forming reaction
                    if not has_alcohol_in_reactants or is_alcohol_forming_reaction:
                        print("Confirmed alcohol formation in final step")
                        late_stage_alcohol = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Late stage alcohol formation detected: {late_stage_alcohol}")
    return late_stage_alcohol

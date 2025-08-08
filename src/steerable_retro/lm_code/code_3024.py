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
    This function detects a strategy where a thiomethyl group is converted to a thione.
    In retrosynthetic analysis, we look for reactions where a thione is converted to a thiomethyl.
    """
    # Track if we found the key feature
    found_thione_formation = False

    def dfs_traverse(node):
        nonlocal found_thione_formation

        if node["type"] == "reaction" and not found_thione_formation:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")

                # In retrosynthesis, we're looking for thione in reactants and thiomethyl in products
                # Check for thiocarbonyl (C=S) in reactants
                has_thiocarbonyl_in_reactants = any(
                    checker.check_fg("Thiocarbonyl", r) for r in reactants_smiles
                )
                print(f"Has thiocarbonyl in reactants: {has_thiocarbonyl_in_reactants}")

                # Check for thioamide in reactants
                has_thioamide_in_reactants = any(
                    checker.check_fg("Thioamide", r) for r in reactants_smiles
                )
                print(f"Has thioamide in reactants: {has_thioamide_in_reactants}")

                # Check for thiourea in reactants
                has_thiourea_in_reactants = any(
                    checker.check_fg("Thiourea", r) for r in reactants_smiles
                )
                print(f"Has thiourea in reactants: {has_thiourea_in_reactants}")

                # Check for thiomethyl group in product
                # Specifically looking for S-CH3, which is a specific type of monosulfide
                has_monosulfide_in_product = checker.check_fg("Monosulfide", product_smiles)
                print(f"Has monosulfide in product: {has_monosulfide_in_product}")

                # Check for relevant reactions that might convert thione to thiomethyl
                relevant_reactions = ["Methylation", "Methylation with MeI_SH", "S-methylation"]

                is_relevant_reaction = any(
                    checker.check_reaction(rxn, rsmi) for rxn in relevant_reactions
                )
                print(f"Is relevant reaction: {is_relevant_reaction}")

                # Consider the strategy found if:
                # 1. Reactants have thiocarbonyl or specific thione-containing groups
                # 2. Product has monosulfide (which could include thiomethyl)
                # 3. The reaction is of a relevant type or we have strong evidence of the conversion
                if (
                    has_thiocarbonyl_in_reactants
                    or has_thioamide_in_reactants
                    or has_thiourea_in_reactants
                ) and has_monosulfide_in_product:
                    print(f"Found thione to thiomethyl conversion in reaction: {rsmi}")
                    found_thione_formation = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_thione_formation

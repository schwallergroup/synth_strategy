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
    This function detects if the synthetic route involves chlorosulfonation of an aromatic ring.
    """
    chlorosulfonation_found = False

    def dfs_traverse(node):
        nonlocal chlorosulfonation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # First check if this is directly classified as aromatic sulfonyl chlorination
            if checker.check_reaction("Aromatic sulfonyl chlorination", rsmi):
                print("Detected aromatic sulfonyl chlorination reaction directly")
                chlorosulfonation_found = True
                return

            # If not directly classified, check for the reaction components
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            # Split reactants and analyze each
            reactants = reactants_str.split(".")

            # Check for chlorosulfonating agents in reactants
            has_chlorosulfonating_agent = False
            has_aromatic_reactant = False
            aromatic_reactant_smiles = ""
            has_sulfonyl_chloride_product = False

            # Common chlorosulfonating agents
            chlorosulfonating_agents = [
                "ClSO3H",
                "O=S(=O)(O)Cl",
                "ClS(=O)(=O)OH",
                "HSO3Cl",  # Chlorosulfonic acid
                "O=S(=O)(Cl)Cl",
                "ClS(=O)(=O)Cl",  # Sulfuryl chloride (SO2Cl2)
                "SOCl2",
                "ClS(=O)Cl",
                "O=S(Cl)Cl",  # Thionyl chloride
            ]

            for reactant in reactants:
                try:
                    # Check if this reactant is a chlorosulfonating agent
                    for agent in chlorosulfonating_agents:
                        if agent in reactant:
                            print(f"Found chlorosulfonating agent: {reactant}")
                            has_chlorosulfonating_agent = True
                            break

                    # Check if this reactant has an aromatic ring
                    if (
                        checker.check_ring("benzene", reactant)
                        or checker.check_ring("naphthalene", reactant)
                        or checker.check_ring("anthracene", reactant)
                    ):
                        print(f"Found aromatic reactant: {reactant}")
                        has_aromatic_reactant = True
                        aromatic_reactant_smiles = reactant
                except Exception as e:
                    print(f"Error processing reactant {reactant}: {e}")

            # Check if product has sulfonyl chloride group
            try:
                if checker.check_fg("Sulfonyl halide", product_str):
                    print(f"Found sulfonyl chloride in product: {product_str}")
                    has_sulfonyl_chloride_product = True
            except Exception as e:
                print(f"Error checking product for sulfonyl chloride: {e}")

            # Pattern-based detection: check if an aromatic ring in reactant gains a sulfonyl chloride in product
            if has_aromatic_reactant and has_sulfonyl_chloride_product:
                # Even if we don't explicitly identify the chlorosulfonating agent, the transformation is clear
                if has_chlorosulfonating_agent:
                    print("Detected chlorosulfonation with identified chlorosulfonating agent")
                    chlorosulfonation_found = True
                else:
                    # Check if the product contains both an aromatic ring and a sulfonyl chloride
                    # This is a strong indicator of chlorosulfonation even if reagent isn't identified
                    print("Detected likely chlorosulfonation based on product formation")
                    chlorosulfonation_found = True

        # Traverse children
        for child in node.get("children", []):
            if not chlorosulfonation_found:  # Stop traversal if already found
                dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return chlorosulfonation_found

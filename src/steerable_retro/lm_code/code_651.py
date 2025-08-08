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
    This function detects the formation of a quinazoline scaffold in the synthesis route.
    """
    quinazoline_formed = False

    def dfs_traverse(node):
        nonlocal quinazoline_formed

        if node["type"] == "reaction" and not quinazoline_formed:
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")

                # Direct check for quinazoline formation
                product_has_quinazoline = checker.check_ring("quinazoline", product_smiles)
                if product_has_quinazoline:
                    print("Product has quinazoline structure")
                    reactants_have_quinazoline = any(
                        checker.check_ring("quinazoline", r) for r in reactants_smiles
                    )

                    if not reactants_have_quinazoline:
                        print("Confirmed: Quinazoline is formed in this reaction")
                        quinazoline_formed = True
                        return

                # Check if this is a Niementowski quinazoline formation reaction
                if checker.check_reaction("Niementowski_quinazoline", rsmi):
                    print("Detected Niementowski_quinazoline reaction")
                    quinazoline_formed = True
                    return

                # Fallback detection method - check for structural changes
                if not quinazoline_formed:
                    # Check if product contains a pyrimidine ring (part of quinazoline)
                    product_has_pyrimidine = checker.check_ring("pyrimidine", product_smiles)
                    product_has_benzene = checker.check_ring("benzene", product_smiles)

                    if product_has_pyrimidine and product_has_benzene:
                        print("Product has pyrimidine and benzene rings (potential quinazoline)")
                        # Check if reactants have pyrimidine
                        reactants_have_pyrimidine = any(
                            checker.check_ring("pyrimidine", r) for r in reactants_smiles
                        )

                        if not reactants_have_pyrimidine:
                            # Additional check for quinoline-like structure (similar to quinazoline)
                            product_has_quinoline = checker.check_ring("quinoline", product_smiles)
                            reactants_have_quinoline = any(
                                checker.check_ring("quinoline", r) for r in reactants_smiles
                            )

                            if product_has_quinoline and not reactants_have_quinoline:
                                print(
                                    "Product has quinoline-like structure not present in reactants"
                                )
                                print(
                                    "Product has both pyrimidine and benzene rings - likely quinazoline"
                                )
                                quinazoline_formed = True

                # Check for pteridin which contains the quinazoline core
                if not quinazoline_formed:
                    product_has_pteridin = checker.check_ring("pteridin", product_smiles)
                    if product_has_pteridin:
                        print("Product has pteridin structure (contains quinazoline core)")
                        reactants_have_pteridin = any(
                            checker.check_ring("pteridin", r) for r in reactants_smiles
                        )

                        if not reactants_have_pteridin:
                            print("Confirmed: Quinazoline core is formed in this reaction")
                            quinazoline_formed = True

            except Exception as e:
                # Handle any errors in extracting or processing reaction data
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            if not quinazoline_formed:  # Stop traversal if quinazoline is already found
                dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return quinazoline_formed

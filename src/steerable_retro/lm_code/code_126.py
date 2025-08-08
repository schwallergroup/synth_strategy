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
    This function detects a synthetic strategy involving the formation of
    a heterocyclic ring system (specifically looking for benzofuran formation).
    """
    # Track if we found heterocycle formation
    found_heterocycle_formation = False

    # List of heterocyclic rings to check for formation
    heterocycles = [
        "benzofuran",  # Primary target
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "indole",
        "furan",
        "pyrrole",
        "thiophene",
        "oxazole",
        "thiazole",
        "imidazole",
    ]

    # List of known heterocycle formation reaction names
    heterocycle_reactions = [
        "benzofuran",
        "{benzofuran}",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "{benzoxazole}",
        "benzothiazole",
        "{benzothiazole}",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "{benzimidazole}",
        "indole",
        "{indole}",
        "Fischer indole",
        "Paal-Knorr pyrrole",
        "{Paal-Knorr pyrrole}",
        "thiazole",
        "{thiazole}",
        "oxazole",
        "imidazole",
        "{imidazole}",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycle_formation

        # If we already found a heterocycle formation, no need to continue
        if found_heterocycle_formation:
            return

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Check if this is a known heterocycle formation reaction
                for reaction_name in heterocycle_reactions:
                    if checker.check_reaction(reaction_name, rsmi):
                        print(f"Found {reaction_name} formation reaction at depth {depth}")
                        found_heterocycle_formation = True
                        return  # Early return once found

                try:
                    # In retrosynthesis, the product is what we're breaking down (left side of SMILES)
                    # and reactants are the precursors (right side of SMILES)
                    product_smiles = rsmi.split(">")[0]
                    reactants_smiles = rsmi.split(">")[-1]
                    reactants_list = reactants_smiles.split(".")

                    # Check for heterocycle formation
                    for heterocycle in heterocycles:
                        # Check if product has the heterocycle
                        product_has_heterocycle = checker.check_ring(heterocycle, product_smiles)

                        if product_has_heterocycle:
                            # Check if any reactant has the heterocycle
                            reactants_have_heterocycle = False
                            for r in reactants_list:
                                if r and checker.check_ring(heterocycle, r):
                                    reactants_have_heterocycle = True
                                    break

                            # In retrosynthesis, if heterocycle is in product (starting material)
                            # but not in reactants (precursors), it means the heterocycle was broken,
                            # which means it was formed in the forward direction
                            if not reactants_have_heterocycle:
                                print(
                                    f"Found heterocycle ({heterocycle}) formation at depth {depth}"
                                )
                                found_heterocycle_formation = True
                                return  # Early return once found
                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Heterocycle formation detected: {found_heterocycle_formation}")
    return found_heterocycle_formation

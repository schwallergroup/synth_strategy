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
    Detects a synthetic strategy involving the coupling of two heterocyclic systems.
    """
    found_heterocycle_coupling = False

    # List of common heterocycles to check
    heterocycles = [
        "furan",
        "pyrrole",
        "thiophene",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "indole",
        "benzofuran",
        "benzothiophene",
        "quinoline",
        "isoquinoline",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
    ]

    # List of common coupling reactions
    coupling_reactions = [
        "Suzuki",
        "Buchwald-Hartwig",
        "N-arylation",
        "Stille",
        "Negishi",
        "Heck",
        "Sonogashira",
        "Ullmann-Goldberg",
        "Ullmann condensation",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycle_coupling

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a coupling reaction
                is_coupling_reaction = False
                for reaction_type in coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_coupling_reaction = True
                        print(f"Found {reaction_type} reaction at depth {depth}")
                        break

                # Check all reactions, but apply different criteria based on depth and reaction type
                # Count heterocycles in reactants
                heterocycle_reactants = []
                for r in reactants:
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, r):
                            heterocycle_reactants.append((r, heterocycle))
                            print(f"Found {heterocycle} in reactant at depth {depth}")
                            break

                # Check if we have at least 2 different heterocycle-containing reactants
                if len(heterocycle_reactants) >= 2:
                    # Verify the product contains heterocycles
                    heterocycles_in_product = []
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, product):
                            heterocycles_in_product.append(heterocycle)
                            print(f"Found {heterocycle} in product at depth {depth}")

                    # If we have heterocycles in reactants and product
                    if len(heterocycles_in_product) >= 1:
                        # For early stages (depth > 1), require a known coupling reaction
                        if depth <= 1 or is_coupling_reaction:
                            found_heterocycle_coupling = True
                            print(f"Found heterocycle coupling at depth {depth}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_heterocycle_coupling

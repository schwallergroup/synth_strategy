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
    This function detects a linear synthesis strategy involving nitrogen heterocycles.

    A linear nitrogen heterocycle strategy involves:
    1. At least 2 reactions
    2. Formation or modification of nitrogen heterocycles
    3. A linear progression of the core N-heterocycle structure
    """
    reaction_count = 0
    n_heterocycle_reactions = 0
    n_heterocycle_rings = [
        "pyrrole",
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
        "pyrrolidine",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "carbazole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "indazole",
        "benzotriazole",
    ]

    # Reactions that form or modify N-heterocycles
    n_heterocycle_forming_reactions = [
        "Paal-Knorr pyrrole synthesis",
        "Fischer indole",
        "Friedlaender chinoline",
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "Niementowski_quinazoline",
        "tetrazole_terminal",
        "tetrazole_connect_regioisomere_1",
        "tetrazole_connect_regioisomere_2",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "3-nitrile-pyridine",
        "pyrazole",
        "indole",
        "imidazole",
    ]

    # Track the main synthetic pathway
    main_pathway = []

    # Track products at each step to establish connectivity
    product_to_reaction = {}

    def dfs_traverse(node, depth=0, path=None):
        nonlocal reaction_count, n_heterocycle_reactions

        if path is None:
            path = []

        current_path = path.copy()

        if node["type"] == "reaction":
            reaction_count += 1

            # Extract product and reactants
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]
            reactants_smiles = rsmi.split(">")[0].split(".")

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if product has nitrogen heterocycle
            product_has_n_heterocycle = False
            product_n_heterocycles = []
            for ring in n_heterocycle_rings:
                if checker.check_ring(ring, product_smiles):
                    product_has_n_heterocycle = True
                    product_n_heterocycles.append(ring)
                    print(f"Product contains {ring}")

            # Check if any reactant has nitrogen heterocycle
            reactants_with_n_heterocycle = 0
            reactant_n_heterocycles = []
            for reactant in reactants_smiles:
                reactant_rings = []
                for ring in n_heterocycle_rings:
                    if checker.check_ring(ring, reactant):
                        reactants_with_n_heterocycle += 1
                        reactant_rings.append(ring)
                        print(f"Reactant contains {ring}: {reactant}")
                if reactant_rings:
                    reactant_n_heterocycles.append((reactant, reactant_rings))

            # Check if this is a N-heterocycle forming or modifying reaction
            is_n_heterocycle_reaction = False

            # Check if reaction is a known N-heterocycle forming reaction
            for rxn_type in n_heterocycle_forming_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    is_n_heterocycle_reaction = True
                    print(f"Detected N-heterocycle reaction: {rxn_type}")
                    break

            # If not a known reaction type, check if it's forming or modifying N-heterocycles
            if not is_n_heterocycle_reaction:
                # Case 1: N-heterocycle is formed (product has it, reactants don't)
                if product_has_n_heterocycle and reactants_with_n_heterocycle == 0:
                    is_n_heterocycle_reaction = True
                    print("N-heterocycle formation detected")

                # Case 2: N-heterocycle is modified (both product and at least one reactant have it)
                elif product_has_n_heterocycle and reactants_with_n_heterocycle > 0:
                    is_n_heterocycle_reaction = True
                    print("N-heterocycle modification detected")

            if is_n_heterocycle_reaction:
                n_heterocycle_reactions += 1
                main_pathway.append((depth, rsmi, product_smiles, product_n_heterocycles))

                # Store the product for connectivity tracking
                product_to_reaction[product_smiles] = (depth, rsmi, product_n_heterocycles)

                # Add to current path
                current_path.append((depth, rsmi, product_smiles, product_n_heterocycles))

        # For molecule nodes, check if they're part of our pathway
        elif node["type"] == "mol" and not node.get("in_stock", False):
            mol_smiles = node["smiles"]
            # Check if this molecule is a product of a reaction in our main pathway
            if mol_smiles in product_to_reaction:
                current_path.append(product_to_reaction[mol_smiles])

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

    # Start traversal from root
    dfs_traverse(route)

    # Sort main pathway by depth
    main_pathway.sort(key=lambda x: x[0])

    print(f"Main pathway reactions: {len(main_pathway)}")
    for item in main_pathway:
        print(f"Depth {item[0]}: {item[1][:50]}... with heterocycles: {item[3]}")

    # Check if we have a connected sequence of N-heterocycle reactions
    is_linear = False
    if len(main_pathway) >= 2:
        # We have at least 2 reactions involving N-heterocycles
        # Check if there's a progression in the heterocycle structure
        # For simplicity, we'll consider it linear if the same types of heterocycles
        # are consistently involved throughout the pathway

        # Get all heterocycle types involved
        all_heterocycles = set()
        for _, _, _, heterocycles in main_pathway:
            all_heterocycles.update(heterocycles)

        # Check if each reaction involves at least one common heterocycle type
        common_heterocycle = False
        for heterocycle in all_heterocycles:
            # Count how many reactions involve this heterocycle
            count = sum(1 for _, _, _, heterocycles in main_pathway if heterocycle in heterocycles)
            if (
                count >= len(main_pathway) * 0.5
            ):  # At least half of the reactions involve this heterocycle
                common_heterocycle = True
                print(
                    f"Common heterocycle found: {heterocycle} in {count}/{len(main_pathway)} reactions"
                )
                break

        is_linear = common_heterocycle

    # Return True if it's a linear synthesis with N heterocycles throughout
    result = reaction_count >= 2 and n_heterocycle_reactions >= 2 and is_linear
    print(
        f"Reaction count: {reaction_count}, N-heterocycle reactions: {n_heterocycle_reactions}, Is linear: {is_linear}"
    )
    return result

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
    This function detects a synthetic strategy involving functionalization of a
    piperazine scaffold at multiple positions.
    """
    # Track if piperazine is present in the final molecule
    piperazine_in_final = False

    # Track functionalization reactions
    functionalization_reactions = []

    def dfs_traverse(node, depth=0):
        nonlocal piperazine_in_final

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for piperazine scaffold in the final molecule (depth 0)
            if depth == 0:
                if checker.check_ring("piperazine", mol_smiles):
                    piperazine_in_final = True
                    print(f"Found piperazine scaffold in final molecule: {mol_smiles}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            try:
                # Extract reactants and product
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check if product contains piperazine
                if checker.check_ring("piperazine", product_part):
                    print(f"Found piperazine in product at depth {depth}: {product_part}")

                    # Check if any reactant contains piperazine
                    piperazine_in_reactants = False
                    for reactant in reactants:
                        if checker.check_ring("piperazine", reactant):
                            piperazine_in_reactants = True
                            break

                    # If piperazine is in both reactants and product, it's likely a functionalization
                    if piperazine_in_reactants:
                        # Check for common functionalization reactions
                        functionalization_types = [
                            "N-alkylation of primary amines with alkyl halides",
                            "N-alkylation of secondary amines with alkyl halides",
                            "Acylation of primary amines",
                            "Acylation of secondary amines",
                            "Reductive amination with aldehyde",
                            "Reductive amination with ketone",
                            "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                            "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
                        ]

                        for rxn_type in functionalization_types:
                            if checker.check_reaction(rxn_type, rsmi):
                                print(f"Found {rxn_type} reaction at depth {depth}")
                                functionalization_reactions.append((depth, rxn_type, rsmi))
                                break

                        # If no specific reaction type was identified, check for general N-functionalization
                        if (
                            len(functionalization_reactions) == 0
                            or functionalization_reactions[-1][0] != depth
                        ):
                            # Check if any nitrogen-containing functional groups were added
                            fg_types = [
                                "Primary amine",
                                "Secondary amine",
                                "Tertiary amine",
                                "Primary amide",
                                "Secondary amide",
                                "Tertiary amide",
                                "Sulfonamide",
                                "Urea",
                                "Carbamate",
                            ]

                            product_fgs = set()
                            for fg in fg_types:
                                if checker.check_fg(fg, product_part):
                                    product_fgs.add(fg)

                            reactant_fgs = set()
                            for reactant in reactants:
                                for fg in fg_types:
                                    if checker.check_fg(fg, reactant):
                                        reactant_fgs.add(fg)

                            # If product has functional groups not in reactants, it's a functionalization
                            if product_fgs - reactant_fgs:
                                print(f"Found general N-functionalization at depth {depth}")
                                functionalization_reactions.append(
                                    (depth, "General N-functionalization", rsmi)
                                )
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Count unique depths where functionalization occurred
    unique_depths = len(set(depth for depth, _, _ in functionalization_reactions))

    print(
        f"Piperazine in final: {piperazine_in_final}, Unique functionalization steps: {unique_depths}"
    )

    # Strategy requires piperazine scaffold with multiple functionalizations (at different depths)
    return piperazine_in_final and unique_depths >= 2

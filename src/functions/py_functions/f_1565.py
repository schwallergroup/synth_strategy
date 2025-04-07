#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects if the synthesis route involves N-methylation pattern
    (presence of N-methyl group in the final product and an N-methylation reaction in the route).
    """
    # Track if we found N-methyl in final product
    n_methyl_in_product = False
    # Track if we found N-methylation reaction
    n_methylation_reaction_found = False

    # List of N-methylation reaction types
    n_methylation_reactions = [
        "N-methylation",
        "Eschweiler-Clarke Primary Amine Methylation",
        "Eschweiler-Clarke Secondary Amine Methylation",
        "Reductive methylation of primary amine with formaldehyde",
        "Methylation with MeI_primary",
        "Methylation with MeI_secondary",
        "Methylation with MeI_tertiary",
        "DMS Amine methylation",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal n_methyl_in_product, n_methylation_reaction_found

        # Check final product for N-methyl groups
        if node["type"] == "mol" and depth == 0:  # Final product
            if "smiles" in node:
                smiles = node["smiles"]
                print(f"Checking final product: {smiles}")

                # Check specifically for N-methyl groups
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Pattern for N-methyl: nitrogen connected to a methyl group
                    pattern = Chem.MolFromSmarts("[N;!$(N=*)]-[CH3]")
                    if mol.HasSubstructMatch(pattern):
                        print("N-methyl group detected in final product")
                        n_methyl_in_product = True

                    # Also check using functional group checker
                    if checker.check_fg("Tertiary amine", smiles):
                        print("Tertiary amine detected in final product")
                        n_methyl_in_product = True

        # Check for N-methylation reactions
        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rxn_smiles = node["metadata"]["rsmi"]
                print(f"Checking reaction: {rxn_smiles}")

                # Check for known N-methylation reactions
                for rxn_type in n_methylation_reactions:
                    if checker.check_reaction(rxn_type, rxn_smiles):
                        print(f"N-methylation reaction detected: {rxn_type}")
                        n_methylation_reaction_found = True
                        break

                # If no specific N-methylation reaction was found, check if the reaction preserves N-methyl groups
                if not n_methylation_reaction_found:
                    try:
                        reactants = rxn_smiles.split(">")[0].split(".")
                        product = rxn_smiles.split(">")[-1]

                        # Check if product has N-methyl
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            pattern = Chem.MolFromSmarts("[N;!$(N=*)]-[CH3]")
                            if product_mol.HasSubstructMatch(pattern):
                                print("Product contains N-methyl group")

                                # Check if any reactant has N-methyl
                                for r in reactants:
                                    if r:
                                        r_mol = Chem.MolFromSmiles(r)
                                        if r_mol and r_mol.HasSubstructMatch(pattern):
                                            print("Reaction preserves N-methyl group")
                                            n_methylation_reaction_found = True
                                            break
                    except Exception as e:
                        print(
                            f"Error analyzing reaction for N-methyl preservation: {e}"
                        )

                # Also check for reactions that convert secondary amines to tertiary amines
                if not n_methylation_reaction_found:
                    try:
                        reactants = rxn_smiles.split(">")[0].split(".")
                        product = rxn_smiles.split(">")[-1]

                        # Check if product has tertiary amine
                        if checker.check_fg("Tertiary amine", product):
                            # Check if any reactant has secondary amine
                            for r in reactants:
                                if r and checker.check_fg("Secondary amine", r):
                                    print(
                                        "Potential N-methylation: Secondary amine to tertiary amine"
                                    )
                                    n_methylation_reaction_found = True
                                    break
                    except Exception as e:
                        print(
                            f"Error analyzing reaction for tertiary amine formation: {e}"
                        )

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both N-methyl in product and N-methylation reaction are found
    result = n_methyl_in_product and n_methylation_reaction_found
    print(
        f"N-methyl in product: {n_methyl_in_product}, N-methylation reaction found: {n_methylation_reaction_found}"
    )
    print(f"Final result: {result}")
    return result

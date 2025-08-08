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
    This function detects a linear synthesis strategy involving an amino acid derivative.
    """
    # Track if amino acid is present in the main synthetic pathway
    amino_acid_in_pathway = False
    # Track if synthesis is linear (no significant branching)
    is_linear = True
    # Track the main synthetic pathway
    main_path_nodes = []

    def dfs_traverse(node, depth=0, path=None):
        nonlocal amino_acid_in_pathway, is_linear, main_path_nodes

        if path is None:
            path = []

        current_path = path + [node]

        # For molecule nodes, check if it's an amino acid
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for amino acid pattern using the checker function
            is_amino_acid = False
            try:
                # Check for carboxylic acid and amine groups in proximity
                has_carboxylic_acid = checker.check_fg("Carboxylic acid", mol_smiles)
                has_amine = (
                    checker.check_fg("Primary amine", mol_smiles)
                    or checker.check_fg("Secondary amine", mol_smiles)
                    or checker.check_fg("Tertiary amine", mol_smiles)
                )

                # More specific check for amino acid structure
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol and has_carboxylic_acid and has_amine:
                    # Check if the amine and carboxylic acid are on adjacent carbons
                    amino_acid_pattern = Chem.MolFromSmarts("[NX3]-[CX4]-[CX3](=[OX1])-[OX2]")
                    if mol.HasSubstructMatch(amino_acid_pattern):
                        is_amino_acid = True
                        print(f"Amino acid detected: {mol_smiles}")
            except Exception as e:
                print(f"Error checking amino acid: {e}")

            # If this is the target molecule (depth 0) or a significant intermediate
            # and it's an amino acid, mark it
            if depth <= 2 and is_amino_acid:
                amino_acid_in_pathway = True
                if node not in main_path_nodes:
                    main_path_nodes.append(node)

        # For reaction nodes, check if it's a typical amino acid modification reaction
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for common amino acid modification reactions
            try:
                is_amino_acid_reaction = (
                    checker.check_reaction("Esterification of Carboxylic Acids", rxn_smiles)
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rxn_smiles
                    )
                    or checker.check_reaction("Boc amine protection", rxn_smiles)
                    or checker.check_reaction("Boc amine deprotection", rxn_smiles)
                    or checker.check_reaction("Reductive amination with aldehyde", rxn_smiles)
                    or checker.check_reaction("Reductive amination with ketone", rxn_smiles)
                )

                if is_amino_acid_reaction:
                    print(f"Amino acid modification reaction detected: {rxn_smiles}")
            except Exception as e:
                print(f"Error checking reaction: {e}")

        children = node.get("children", [])

        # If this node has more than 2 children, it might indicate a non-linear synthesis
        # But we need to check if these are just reagents or actual branching pathways
        if len(children) > 2:
            significant_children = 0
            for child in children:
                if child["type"] == "mol" and not child.get("in_stock", False):
                    significant_children += 1

            if significant_children > 1:
                print(
                    f"Non-linear synthesis detected: node has {significant_children} significant children"
                )
                is_linear = False

        # Continue traversal
        for child in children:
            dfs_traverse(child, depth + 1, current_path)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have a linear amino acid synthesis
    result = amino_acid_in_pathway and is_linear

    if result:
        print("Linear amino acid synthesis strategy detected")
    else:
        if not amino_acid_in_pathway:
            print("No amino acid detected in the main synthetic pathway")
        if not is_linear:
            print("Synthesis is not linear")

    return result

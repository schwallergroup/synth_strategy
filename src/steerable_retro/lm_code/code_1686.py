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


def main(route):
    """
    This function detects if the synthetic route includes a deprotection step
    (acetamide to amine) followed by an amidation reaction.
    """
    # Track reactions in sequence
    reaction_sequence = []

    def dfs_traverse(node):
        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check reaction type
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol:
                # Check for deprotection (acetamide to amine)
                acetamide_pattern = Chem.MolFromSmarts("[#6][#6](=[#8])[#7]")
                amine_pattern = Chem.MolFromSmarts("[#7;H2]")

                is_deprotection = False
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if (
                        reactant_mol
                        and reactant_mol.HasSubstructMatch(acetamide_pattern)
                        and product_mol.HasSubstructMatch(amine_pattern)
                    ):
                        is_deprotection = True
                        reaction_sequence.append("deprotection")
                        print("Detected acetamide deprotection")
                        break

                # Check for amidation
                if not is_deprotection:
                    amide_pattern = Chem.MolFromSmarts("[#6](=[#8])[#7]")
                    if product_mol.HasSubstructMatch(amide_pattern):
                        # Check if this is a new amide formation
                        amide_count_product = len(product_mol.GetSubstructMatches(amide_pattern))
                        amide_count_reactants = 0

                        for reactant in reactants_smiles:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                amide_count_reactants += len(
                                    reactant_mol.GetSubstructMatches(amide_pattern)
                                )

                        if amide_count_product > amide_count_reactants:
                            reaction_sequence.append("amidation")
                            print("Detected amide formation")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if deprotection is followed by amidation
    for i in range(len(reaction_sequence) - 1):
        if reaction_sequence[i] == "deprotection" and reaction_sequence[i + 1] == "amidation":
            return True

    return False

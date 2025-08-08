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


def main(route):
    """
    Detects if the synthesis involves multiple aromatic substitution reactions
    (bromination, nitration, etc.) on the same core structure.
    """
    aromatic_substitutions = 0

    def dfs_traverse(node):
        nonlocal aromatic_substitutions

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Patterns for common aromatic substitutions
                bromination_pattern = Chem.MolFromSmarts("c-Br")
                nitration_pattern = Chem.MolFromSmarts("c-[N+](=O)[O-]")

                product_mol = Chem.MolFromSmiles(product_smiles)
                reactant_mols = [
                    Chem.MolFromSmiles(r) for r in reactants_smiles if Chem.MolFromSmiles(r)
                ]

                # Check if product has a new aromatic substitution not present in reactants
                if product_mol:
                    # Check for new bromination
                    if product_mol.HasSubstructMatch(bromination_pattern):
                        br_in_product = len(product_mol.GetSubstructMatches(bromination_pattern))
                        br_in_reactants = sum(
                            [
                                len(r.GetSubstructMatches(bromination_pattern))
                                for r in reactant_mols
                                if r
                            ]
                        )

                        if br_in_product > br_in_reactants:
                            print("Aromatic bromination detected")
                            aromatic_substitutions += 1

                    # Check for new nitration
                    if product_mol.HasSubstructMatch(nitration_pattern):
                        nitro_in_product = len(product_mol.GetSubstructMatches(nitration_pattern))
                        nitro_in_reactants = sum(
                            [
                                len(r.GetSubstructMatches(nitration_pattern))
                                for r in reactant_mols
                                if r
                            ]
                        )

                        if nitro_in_product > nitro_in_reactants:
                            print("Aromatic nitration detected")
                            aromatic_substitutions += 1
            except Exception as e:
                print(f"Error in multiple_aromatic_substitution_strategy: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least 2 aromatic substitutions are detected
    return aromatic_substitutions >= 2

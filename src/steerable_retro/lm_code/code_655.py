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
    This function detects if a bromo substituent is preserved throughout the synthesis.
    It checks if the target molecule has a bromine and if this bromine is preserved
    in all non-stock intermediates throughout the synthesis route.
    """
    # First check if the target molecule has a bromine
    target_mol_smiles = route["smiles"]
    target_mol = Chem.MolFromSmiles(target_mol_smiles)
    bromo_pattern = Chem.MolFromSmarts("[Br]")

    if not target_mol.HasSubstructMatch(bromo_pattern):
        print(f"Target molecule doesn't have bromine: {target_mol_smiles}")
        return False

    # Track bromine preservation through the synthesis
    # In retrosynthesis, we need to track when bromine is introduced
    # and ensure it's preserved from that point forward

    # Store nodes where bromine should be present
    nodes_with_bromine = set()

    def dfs_first_pass(node, path=None):
        """First pass to identify where bromine is introduced and should be preserved"""
        if path is None:
            path = []

        current_path = path + [node]

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            mol = Chem.MolFromSmiles(mol_smiles)

            # If this molecule has bromine, all molecules in the path to the target should have bromine
            if mol and mol.HasSubstructMatch(bromo_pattern):
                for n in current_path:
                    if n["type"] == "mol" and not n.get("in_stock", False):
                        nodes_with_bromine.add(id(n))

        # Traverse children
        for child in node.get("children", []):
            dfs_first_pass(child, current_path)

    # Second pass to check if bromine is preserved where it should be
    def check_bromine_preservation():
        bromine_preserved = True

        def dfs_check(node):
            nonlocal bromine_preserved

            if node["type"] == "mol":
                # If this node should have bromine (based on first pass)
                if id(node) in nodes_with_bromine:
                    mol_smiles = node["smiles"]
                    mol = Chem.MolFromSmiles(mol_smiles)
                    if mol and not mol.HasSubstructMatch(bromo_pattern):
                        print(f"Intermediate without bromine detected: {mol_smiles}")
                        bromine_preserved = False

            elif node["type"] == "reaction":
                try:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    product_mol = Chem.MolFromSmiles(product)
                    product_has_bromo = product_mol and product_mol.HasSubstructMatch(bromo_pattern)

                    # Check if any reactant has bromine
                    reactant_has_bromo = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(bromo_pattern):
                            reactant_has_bromo = True
                            break

                    # In retrosynthesis, if product has bromine but reactants don't,
                    # it means bromine is introduced at this step in forward synthesis
                    if product_has_bromo and not reactant_has_bromo:
                        print(f"Bromine is introduced in this reaction: {rsmi}")

                    # If product has bromine but reactants don't have it in a reaction where
                    # bromine should be preserved, that's an issue
                    if product_has_bromo and not reactant_has_bromo:
                        for child in node.get("children", []):
                            if child["type"] == "mol" and id(child) in nodes_with_bromine:
                                print(f"Bromine preservation issue in reaction: {rsmi}")
                                bromine_preserved = False
                                break

                except Exception as e:
                    print(f"Error analyzing reaction: {e}")

            # Traverse children
            for child in node.get("children", []):
                dfs_check(child)

        dfs_check(route)
        return bromine_preserved

    # Run first pass to identify where bromine should be present
    dfs_first_pass(route)

    # Run second pass to check if bromine is preserved where it should be
    return check_bromine_preservation()

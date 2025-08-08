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
    This function detects if the synthetic route involves nitrile-mediated heterocycle formation,
    where a nitrile group is maintained through multiple steps before being used in cyclization.
    """
    nitrile_present_steps = 0
    nitrile_used_in_cyclization = False

    def dfs_traverse(node):
        nonlocal nitrile_present_steps, nitrile_used_in_cyclization

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(r for r in reactant_mols if r):
                # Check for nitrile presence
                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                pyrazole_pattern = Chem.MolFromSmarts("c1[n]n[c]c1")

                reactants_have_nitrile = any(
                    r.HasSubstructMatch(nitrile_pattern) for r in reactant_mols if r
                )
                product_has_nitrile = product_mol.HasSubstructMatch(nitrile_pattern)

                # Count steps where nitrile is present
                if reactants_have_nitrile or product_has_nitrile:
                    nitrile_present_steps += 1

                # Check if nitrile is used in heterocycle formation
                if (
                    reactants_have_nitrile
                    and not product_has_nitrile
                    and product_mol.HasSubstructMatch(pyrazole_pattern)
                ):
                    print("Nitrile used in heterocycle formation")
                    nitrile_used_in_cyclization = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if nitrile is present in multiple steps and used in cyclization
    return nitrile_present_steps >= 2 and nitrile_used_in_cyclization

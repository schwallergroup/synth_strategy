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
    This function detects if the synthesis uses a convergent approach with N-alkylation
    to join fragments.
    """
    has_convergent_n_alkylation = False

    def dfs_traverse(node):
        nonlocal has_convergent_n_alkylation

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if we have multiple reactants (convergent)
            if len(reactants) >= 2:
                # Check for N-alkylation pattern
                # Secondary amine in one fragment
                sec_amine_pattern = Chem.MolFromSmarts("[NH][#6]")

                # Leaving group in another fragment
                leaving_group_pattern = Chem.MolFromSmarts("[#6][S](=O)(=O)[O][#6]")  # Mesylate

                # Product should have tertiary amine
                tert_amine_pattern = Chem.MolFromSmarts("[N]([#6])([#6])[#6]")

                has_sec_amine = False
                has_leaving_group = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(sec_amine_pattern):
                        has_sec_amine = True
                    if mol and mol.HasSubstructMatch(leaving_group_pattern):
                        has_leaving_group = True

                product_mol = Chem.MolFromSmiles(product)
                has_tert_amine = product_mol and product_mol.HasSubstructMatch(tert_amine_pattern)

                if has_sec_amine and has_leaving_group and has_tert_amine:
                    has_convergent_n_alkylation = True
                    print("Detected convergent N-alkylation strategy")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_convergent_n_alkylation

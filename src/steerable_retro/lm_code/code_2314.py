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
    This function detects if the synthesis involves N-alkylation of a piperazine scaffold.
    """
    n_alkylation_detected = False

    # SMARTS pattern for piperazine
    piperazine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6][#6]1")

    def dfs_traverse(node):
        nonlocal n_alkylation_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")
            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product_part)

            if product_mol and product_mol.HasSubstructMatch(piperazine_pattern):
                # Check for alkyl halide in reactants
                alkyl_halide_pattern = Chem.MolFromSmarts("[#6][Cl,Br,I]")
                piperazine_in_reactants = False
                alkyl_halide_in_reactants = False

                for r_mol in reactants_mols:
                    if r_mol:
                        if r_mol.HasSubstructMatch(piperazine_pattern):
                            piperazine_in_reactants = True
                        if r_mol.HasSubstructMatch(alkyl_halide_pattern):
                            alkyl_halide_in_reactants = True

                if piperazine_in_reactants and alkyl_halide_in_reactants:
                    n_alkylation_detected = True
                    print(f"N-alkylation of piperazine detected in: {rsmi}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    print(f"N-alkylation piperazine strategy: {n_alkylation_detected}")
    return n_alkylation_detected

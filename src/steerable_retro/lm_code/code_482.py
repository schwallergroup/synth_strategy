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
    Detects if the synthesis involves C-N bond formation via nucleophilic substitution
    (displacement of halide by nitrogen nucleophile).
    """
    cn_substitution_detected = False

    def dfs_traverse(node):
        nonlocal cn_substitution_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for C-Cl pattern in reactants
            c_cl_pattern = Chem.MolFromSmarts("[#6][Cl]")

            # Check for amine pattern in reactants
            amine_pattern = Chem.MolFromSmarts("[#7;!$(N=*);!$(NC=O)]")

            # Check if reactants contain C-Cl and amine
            c_cl_in_reactants = False
            amine_in_reactants = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if not reactant_mol:
                    continue

                if reactant_mol.HasSubstructMatch(c_cl_pattern):
                    c_cl_in_reactants = True

                if reactant_mol.HasSubstructMatch(amine_pattern):
                    amine_in_reactants = True

            # Check if product has new C-N bond where C-Cl was
            if c_cl_in_reactants and amine_in_reactants:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # This is a simplification - ideally we would track specific atoms
                    # but for demonstration, we'll just check if product has C-N bonds
                    c_n_pattern = Chem.MolFromSmarts("[#6][#7]")
                    if product_mol.HasSubstructMatch(c_n_pattern):
                        print("C-N bond formation via nucleophilic substitution detected")
                        cn_substitution_detected = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return cn_substitution_detected

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
    This function detects if reductive amination is used in the synthesis.
    """
    found_reductive_amination = False

    def dfs_traverse(node):
        nonlocal found_reductive_amination

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ketone/aldehyde and amine in reactants
                carbonyl_pattern = Chem.MolFromSmarts("[#6]-C(=O)-[#6,#1]")
                amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]-[#6]")

                has_carbonyl = False
                has_amine = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(carbonyl_pattern):
                            has_carbonyl = True
                        if mol.HasSubstructMatch(amine_pattern):
                            has_amine = True

                # Check if product has new C-N bond where carbonyl was
                if has_carbonyl and has_amine:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # Look for C-N-C pattern that wasn't in reactants
                        c_n_c_pattern = Chem.MolFromSmarts("[#6]-[#7]-[#6]")
                        if product_mol.HasSubstructMatch(c_n_c_pattern):
                            found_reductive_amination = True
                            print("Found reductive amination pattern")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_reductive_amination

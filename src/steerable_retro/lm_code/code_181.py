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
    This function detects the installation of a long alkoxy chain onto an aromatic ring.
    """
    found_alkoxy_installation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_alkoxy_installation

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for phenol pattern in reactants
                phenol_pattern = Chem.MolFromSmarts("[c][OH]")

                # Check for long alkyl chain pattern
                long_chain_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2]")

                # Check for alkoxy product pattern
                alkoxy_pattern = Chem.MolFromSmarts("[c][O][CH2][CH2][CH2][CH2][CH2]")

                # Check reactants and product
                has_phenol = False
                has_long_chain = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(phenol_pattern):
                            has_phenol = True
                        if reactant_mol.HasSubstructMatch(long_chain_pattern):
                            has_long_chain = True

                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(alkoxy_pattern):
                    if has_phenol and has_long_chain:
                        found_alkoxy_installation = True
                        print("Found long alkoxy chain installation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_alkoxy_installation

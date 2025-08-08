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
    This function detects a nitro group reduction to amine in the synthesis.
    """
    found_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal found_nitro_reduction

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro reduction
                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol:
                    # Check if product contains amine
                    amine_pattern = Chem.MolFromSmarts("[#7;H2]")

                    # Check if any reactant contains nitro
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")

                    reactants_have_nitro = False
                    for mol in reactant_mols:
                        if mol and mol.HasSubstructMatch(nitro_pattern):
                            reactants_have_nitro = True
                            break

                    if product_mol.HasSubstructMatch(amine_pattern) and reactants_have_nitro:
                        found_nitro_reduction = True
                        print("Found nitro reduction to amine")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_nitro_reduction

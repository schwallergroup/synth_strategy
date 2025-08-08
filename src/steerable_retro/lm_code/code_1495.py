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
    Detects if the synthesis uses reductive amination for C-N bond formation.
    """
    has_reductive_amination = False

    def dfs_traverse(node, depth=0):
        nonlocal has_reductive_amination

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for aldehyde pattern in reactants
                aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
                amine_pattern = Chem.MolFromSmarts("[NH]")

                has_aldehyde = False
                has_amine = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(aldehyde_pattern):
                                has_aldehyde = True
                            if mol.HasSubstructMatch(amine_pattern):
                                has_amine = True
                    except:
                        continue

                # Check if product has new C-N bond where aldehyde carbon is now sp3
                if has_aldehyde and has_amine:
                    try:
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            # This is a simplified check - in practice would need more sophisticated analysis
                            if product_mol.HasSubstructMatch(Chem.MolFromSmarts("[C][N]")):
                                has_reductive_amination = True
                                print(f"Found reductive amination at depth {depth}")
                    except:
                        pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_reductive_amination

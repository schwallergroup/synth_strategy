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
    This function detects diaryl ether formation via SNAr reaction.
    Looks for a reaction where a C-O bond is formed between two aromatic rings.
    """
    diaryl_ether_formed = False

    def dfs_traverse(node):
        nonlocal diaryl_ether_formed

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains diaryl ether
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    diaryl_ether_pattern = Chem.MolFromSmarts("[c]-[O]-[c]")
                    if product_mol.HasSubstructMatch(diaryl_ether_pattern):
                        # Check if reactants include aromatic hydroxyl and fluoroaromatic
                        has_aromatic_oh = False
                        has_fluoroaromatic = False

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[OH]")):
                                    has_aromatic_oh = True
                                if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[F]")):
                                    has_fluoroaromatic = True

                        if has_aromatic_oh and has_fluoroaromatic:
                            diaryl_ether_formed = True
                            print("Detected diaryl ether formation via SNAr")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return diaryl_ether_formed

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
    This function detects if the synthesis involves formation of a sulfonamide
    bond (S-N) from a sulfonyl chloride and an amine.
    """
    has_sulfonamide_formation = False

    def dfs_traverse(node):
        nonlocal has_sulfonamide_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains sulfonamide
                sulfonamide_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])-[#7]")
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(sulfonamide_pattern):
                    # Check if reactants contain sulfonyl chloride and amine
                    sulfonyl_chloride_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])-[#17]")
                    amine_pattern = Chem.MolFromSmarts("[#7;H2]")

                    has_sulfonyl_chloride = False
                    has_amine = False

                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            if mol.HasSubstructMatch(sulfonyl_chloride_pattern):
                                has_sulfonyl_chloride = True
                                print("Found sulfonyl chloride in reactants")
                            if mol.HasSubstructMatch(amine_pattern):
                                has_amine = True
                                print("Found amine in reactants")

                    if has_sulfonyl_chloride and has_amine:
                        has_sulfonamide_formation = True
                        print("Detected sulfonamide formation")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_sulfonamide_formation

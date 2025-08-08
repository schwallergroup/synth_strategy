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
    Detects if the synthetic route includes amide formation from an acid chloride and amine.
    """
    has_amide_formation = False

    def dfs_traverse(node):
        nonlocal has_amide_formation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Split into reactants and product
                parts = rsmi.split(">")
                if len(parts) >= 3:
                    reactants = parts[0].split(".")
                    product = parts[-1]

                    # Check for acid chloride in reactants
                    acid_chloride_pattern = Chem.MolFromSmarts("C(=O)Cl")

                    # Check for amine in reactants
                    amine_pattern = Chem.MolFromSmarts("[NH2]")

                    has_acid_chloride = False
                    has_amine = False

                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                if mol.HasSubstructMatch(acid_chloride_pattern):
                                    has_acid_chloride = True
                                if mol.HasSubstructMatch(amine_pattern):
                                    has_amine = True
                        except:
                            continue

                    # Check if product has an amide bond
                    if has_acid_chloride and has_amine:
                        try:
                            product_mol = Chem.MolFromSmiles(product)
                            amide_pattern = Chem.MolFromSmarts("C(=O)N")
                            if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                                print("Found amide formation from acid chloride and amine")
                                has_amide_formation = True
                        except:
                            pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_amide_formation

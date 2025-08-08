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
    This function detects a synthetic strategy involving amide formation from
    an acid chloride and an amine.
    """
    has_amide_formation = False

    def dfs_traverse(node):
        nonlocal has_amide_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)

                # Check for acid chloride in reactants
                has_acid_chloride = False
                has_amine = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        acid_chloride_pattern = Chem.MolFromSmarts("[C](=[O])[Cl]")
                        amine_pattern = Chem.MolFromSmarts("[NH2][c]")

                        if reactant_mol.HasSubstructMatch(acid_chloride_pattern):
                            has_acid_chloride = True
                        if reactant_mol.HasSubstructMatch(amine_pattern):
                            has_amine = True

                # Check if product has an amide
                if product_mol:
                    amide_pattern = Chem.MolFromSmarts("[C](=[O])[NH][c]")
                    if (
                        product_mol.HasSubstructMatch(amide_pattern)
                        and has_acid_chloride
                        and has_amine
                    ):
                        has_amide_formation = True
                        print(f"Detected amide formation in reaction: {rsmi}")
            except:
                pass

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Amide formation strategy detected: {has_amide_formation}")
    return has_amide_formation

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
    This function detects if the synthetic route contains a sequence of
    functional group interconversions: ester → acid → alcohol → chloride.
    """
    # Track if we've seen each transformation
    ester_to_acid = False
    acid_to_alcohol = False
    alcohol_to_chloride = False

    def dfs_traverse(node):
        nonlocal ester_to_acid, acid_to_alcohol, alcohol_to_chloride

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Patterns for functional groups
            ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")
            acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")
            alcohol_pattern = Chem.MolFromSmarts("[c][CH2][OH]")
            chloromethyl_pattern = Chem.MolFromSmarts("[c][CH2][Cl]")

            try:
                product_mol = Chem.MolFromSmiles(product)

                # Check for ester to acid conversion
                if not ester_to_acid:
                    reactant_has_ester = False
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(ester_pattern):
                                reactant_has_ester = True
                                break
                        except:
                            continue

                    if (
                        reactant_has_ester
                        and product_mol
                        and product_mol.HasSubstructMatch(acid_pattern)
                    ):
                        print("Ester to acid conversion detected")
                        ester_to_acid = True

                # Check for acid to alcohol conversion
                if not acid_to_alcohol:
                    reactant_has_acid = False
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(acid_pattern):
                                reactant_has_acid = True
                                break
                        except:
                            continue

                    if (
                        reactant_has_acid
                        and product_mol
                        and product_mol.HasSubstructMatch(alcohol_pattern)
                    ):
                        print("Acid to alcohol conversion detected")
                        acid_to_alcohol = True

                # Check for alcohol to chloride conversion
                if not alcohol_to_chloride:
                    reactant_has_alcohol = False
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.HasSubstructMatch(alcohol_pattern):
                                reactant_has_alcohol = True
                                break
                        except:
                            continue

                    if (
                        reactant_has_alcohol
                        and product_mol
                        and product_mol.HasSubstructMatch(chloromethyl_pattern)
                    ):
                        print("Alcohol to chloride conversion detected")
                        alcohol_to_chloride = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    # Return True if all three transformations are found
    return ester_to_acid and acid_to_alcohol and alcohol_to_chloride

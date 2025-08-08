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
    Detects if the synthesis involves a specific sequence of functional group transformations:
    carboxylic acid → Weinreb amide → ketone → bromomethyl ketone.
    """
    # Track if we've seen each transformation in the sequence
    seen_transformations = {
        "acid_to_weinreb": False,
        "weinreb_to_ketone": False,
        "ketone_to_bromomethyl": False,
    }

    def dfs_traverse(node):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Define patterns for functional groups
                carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")
                weinreb_amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]([CH3])[O][CH3]")
                ketone_pattern = Chem.MolFromSmarts("[C](=[O])[CH3]")
                bromomethyl_ketone_pattern = Chem.MolFromSmarts("[C](=[O])[CH2][Br]")

                # Check for acid to Weinreb amide transformation
                if not seen_transformations["acid_to_weinreb"]:
                    has_acid = False
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(carboxylic_acid_pattern):
                            has_acid = True
                            break

                    product_mol = Chem.MolFromSmiles(product)
                    if (
                        has_acid
                        and product_mol
                        and product_mol.HasSubstructMatch(weinreb_amide_pattern)
                    ):
                        seen_transformations["acid_to_weinreb"] = True
                        print("Detected carboxylic acid to Weinreb amide transformation")

                # Check for Weinreb amide to ketone transformation
                if not seen_transformations["weinreb_to_ketone"]:
                    has_weinreb = False
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(weinreb_amide_pattern):
                            has_weinreb = True
                            break

                    product_mol = Chem.MolFromSmiles(product)
                    if (
                        has_weinreb
                        and product_mol
                        and product_mol.HasSubstructMatch(ketone_pattern)
                    ):
                        seen_transformations["weinreb_to_ketone"] = True
                        print("Detected Weinreb amide to ketone transformation")

                # Check for ketone to bromomethyl ketone transformation
                if not seen_transformations["ketone_to_bromomethyl"]:
                    has_ketone = False
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(ketone_pattern):
                            has_ketone = True
                            break

                    product_mol = Chem.MolFromSmiles(product)
                    if (
                        has_ketone
                        and product_mol
                        and product_mol.HasSubstructMatch(bromomethyl_ketone_pattern)
                    ):
                        seen_transformations["ketone_to_bromomethyl"] = True
                        print("Detected ketone to bromomethyl ketone transformation")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if we've seen the complete sequence
    return (
        seen_transformations["acid_to_weinreb"]
        and seen_transformations["weinreb_to_ketone"]
        and seen_transformations["ketone_to_bromomethyl"]
    )

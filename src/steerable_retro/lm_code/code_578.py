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
    aldehyde → alkene → alcohol → ether → acid → amide
    """
    # Initialize tracking for each transformation
    transformations = {
        "aldehyde_to_alkene": False,
        "alkene_to_alcohol": False,
        "alcohol_to_ether": False,
        "ester_to_acid": False,
        "acid_to_amide": False,
    }

    def dfs_traverse(node):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Try to create RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]
                product_mol = Chem.MolFromSmiles(product)

                if not product_mol or not reactant_mols:
                    return

                # Check for aldehyde to alkene (Wittig)
                aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")
                alkene_pattern = Chem.MolFromSmarts("C=C")
                if any(
                    mol.HasSubstructMatch(aldehyde_pattern) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(alkene_pattern):
                    transformations["aldehyde_to_alkene"] = True
                    print("Detected aldehyde to alkene transformation")

                # Check for alkene to alcohol
                if any(
                    mol.HasSubstructMatch(alkene_pattern) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]")):
                    transformations["alkene_to_alcohol"] = True
                    print("Detected alkene to alcohol transformation")

                # Check for alcohol to ether
                alcohol_pattern = Chem.MolFromSmarts("[OH]")
                ether_pattern = Chem.MolFromSmarts("[#6]-[O]-[#6]")
                if any(
                    mol.HasSubstructMatch(alcohol_pattern) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(ether_pattern):
                    transformations["alcohol_to_ether"] = True
                    print("Detected alcohol to ether transformation")

                # Check for ester to acid
                ester_pattern = Chem.MolFromSmarts("C(=O)OC")
                acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
                if any(
                    mol.HasSubstructMatch(ester_pattern) for mol in reactant_mols
                ) and product_mol.HasSubstructMatch(acid_pattern):
                    transformations["ester_to_acid"] = True
                    print("Detected ester to acid transformation")

                # Check for acid to amide
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                amide_pattern = Chem.MolFromSmarts("C(=O)[NH]")
                if (
                    any(mol.HasSubstructMatch(acid_pattern) for mol in reactant_mols)
                    and any(mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols)
                    and product_mol.HasSubstructMatch(amide_pattern)
                ):
                    transformations["acid_to_amide"] = True
                    print("Detected acid to amide transformation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if at least 3 of the transformations are present
    transformation_count = sum(1 for value in transformations.values() if value)
    if transformation_count >= 3:
        print(f"Detected {transformation_count} of 5 expected functional group transformations")
        return True
    return False

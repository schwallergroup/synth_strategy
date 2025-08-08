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
    This function detects if the synthetic route involves a sequential
    reduction of carbonyl groups: ester → aldehyde → alcohol.
    """
    # Track if we've seen each transformation
    ester_to_aldehyde = False
    aldehyde_to_alcohol = False

    def dfs_traverse(node):
        nonlocal ester_to_aldehyde, aldehyde_to_alcohol

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if product_mol and reactants_mols:
                    # Check for ester to aldehyde conversion
                    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
                    aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")

                    # Check if any reactant has ester and product has aldehyde
                    for reactant in reactants_mols:
                        if (
                            reactant
                            and reactant.HasSubstructMatch(ester_pattern)
                            and product_mol.HasSubstructMatch(aldehyde_pattern)
                        ):
                            print("Detected ester to aldehyde conversion")
                            ester_to_aldehyde = True
                            break

                    # Check for aldehyde to alcohol conversion
                    alcohol_pattern = Chem.MolFromSmarts("[CH][OH]")

                    # Check if any reactant has aldehyde and product has alcohol
                    for reactant in reactants_mols:
                        if (
                            reactant
                            and reactant.HasSubstructMatch(aldehyde_pattern)
                            and product_mol.HasSubstructMatch(alcohol_pattern)
                        ):
                            print("Detected aldehyde to alcohol conversion")
                            aldehyde_to_alcohol = True
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True only if both transformations are found
    return ester_to_aldehyde and aldehyde_to_alcohol

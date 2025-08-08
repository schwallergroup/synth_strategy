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
    This function detects if the synthetic route involves a specific sequence of
    functional group interconversions: methyl → bromomethyl → methoxymethyl
    """
    # Track the transformations we find
    methyl_to_bromomethyl_found = False
    bromomethyl_to_methoxymethyl_found = False

    def dfs_traverse(node):
        nonlocal methyl_to_bromomethyl_found, bromomethyl_to_methoxymethyl_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(r for r in reactant_mols):
                # Check for methyl to bromomethyl transformation
                methyl_pattern = Chem.MolFromSmarts("[CH3]")
                bromomethyl_pattern = Chem.MolFromSmarts("[CH2][Br]")

                reactants_have_methyl = any(
                    r.HasSubstructMatch(methyl_pattern) for r in reactant_mols if r
                )
                product_has_bromomethyl = product_mol.HasSubstructMatch(bromomethyl_pattern)

                if reactants_have_methyl and product_has_bromomethyl:
                    print("Methyl to bromomethyl transformation detected")
                    methyl_to_bromomethyl_found = True

                # Check for bromomethyl to methoxymethyl transformation
                methoxymethyl_pattern = Chem.MolFromSmarts("[CH2][O][CH3]")

                reactants_have_bromomethyl = any(
                    r.HasSubstructMatch(bromomethyl_pattern) for r in reactant_mols if r
                )
                product_has_methoxymethyl = product_mol.HasSubstructMatch(methoxymethyl_pattern)

                if reactants_have_bromomethyl and product_has_methoxymethyl:
                    print("Bromomethyl to methoxymethyl transformation detected")
                    bromomethyl_to_methoxymethyl_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both transformations were found
    return methyl_to_bromomethyl_found and bromomethyl_to_methoxymethyl_found

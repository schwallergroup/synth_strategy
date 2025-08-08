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
    This function detects if multiple SNAr reactions (nucleophilic aromatic substitution)
    are used in the synthesis, identified by Cl leaving group and N nucleophile.
    """
    snar_count = 0

    def dfs_traverse(node):
        nonlocal snar_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Look for patterns indicating SNAr:
            # 1. Aromatic chloride in one reactant
            # 2. Amine nucleophile in another reactant
            # 3. C-N bond formation in product

            chloro_aromatic_pattern = Chem.MolFromSmarts("c[Cl]")
            amine_pattern = Chem.MolFromSmarts("[N;!$(N=*);!$(N#*);!$(N[N,O,S])]")

            has_chloro_aromatic = False
            has_amine = False

            for reactant_smiles in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                if not reactant_mol:
                    continue

                if reactant_mol.HasSubstructMatch(chloro_aromatic_pattern):
                    has_chloro_aromatic = True

                if reactant_mol.HasSubstructMatch(amine_pattern):
                    has_amine = True

            # If both patterns are found in reactants, check if product has new C-N bond
            if has_chloro_aromatic and has_amine:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    # This is a simplification - a more robust implementation would
                    # check for the specific C-N bond formation at the position where Cl was
                    snar_count += 1
                    print(f"Potential SNAr reaction detected: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Check if multiple SNAr reactions were detected
    multiple_snar = snar_count >= 2

    if multiple_snar:
        print(f"Multiple SNAr reactions detected: {snar_count}")

    return multiple_snar

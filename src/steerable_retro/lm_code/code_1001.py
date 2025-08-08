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
    Detects if the synthetic route uses N-alkylation strategy with chloro-containing intermediates.
    """
    n_alkylation_with_chloro = False
    secondary_amine_pattern = Chem.MolFromSmarts("[NH]")
    tertiary_amine_pattern = Chem.MolFromSmarts("[#7](-[#6])(-[#6])-[#6]")
    chloro_pattern = Chem.MolFromSmarts("[#6]-[Cl]")

    def dfs_traverse(node):
        nonlocal n_alkylation_with_chloro

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            products_smiles = rsmi.split(">")[-1]

            try:
                reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                products_mol = Chem.MolFromSmiles(products_smiles)

                if reactants_mol and products_mol:
                    # Check for secondary amine in reactants, tertiary amine in products, and chloro group in reactants
                    if (
                        reactants_mol.HasSubstructMatch(secondary_amine_pattern)
                        and products_mol.HasSubstructMatch(tertiary_amine_pattern)
                        and reactants_mol.HasSubstructMatch(chloro_pattern)
                    ):
                        n_alkylation_with_chloro = True
                        print(f"Found N-alkylation with chloro intermediate: {rsmi}")
            except:
                print(f"Error processing reaction SMILES: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"N-alkylation strategy with chloro intermediates detected: {n_alkylation_with_chloro}")
    return n_alkylation_with_chloro

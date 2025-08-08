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
    This function detects if the synthetic route involves amide bond formation.
    """
    amide_formation_found = False

    def dfs_traverse(node):
        nonlocal amide_formation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactants contain a carboxylic acid and an amine, and product contains an amide
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and len(reactant_mols) >= 2:
                # Check for carboxylic acid in reactants
                acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H]")
                has_acid = any(mol.HasSubstructMatch(acid_pattern) for mol in reactant_mols if mol)

                # Check for amine in reactants
                amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1]")
                has_amine = any(
                    mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols if mol
                )

                # Check for amide in product
                amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3;H1,H0]")
                has_amide = product_mol.HasSubstructMatch(amide_pattern) if product_mol else False

                if has_acid and has_amine and has_amide:
                    print("Amide formation detected")
                    amide_formation_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return amide_formation_found

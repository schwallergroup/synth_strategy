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
    Detects if the synthesis involves O-alkylation of a phenol to form an ether.
    """
    has_o_alkylation = False

    def dfs_traverse(node):
        nonlocal has_o_alkylation

        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Create molecules
            product_mol = Chem.MolFromSmiles(product)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

            if product_mol and all(reactant_mols):
                # Check for O-alkylation pattern
                phenol_pattern = Chem.MolFromSmarts("c[OH]")
                alcohol_pattern = Chem.MolFromSmarts("[CH]([CH3])[OH]")  # isopropanol pattern
                ether_pattern = Chem.MolFromSmarts("cO[CH]([CH3])[CH3]")  # isopropoxy pattern

                has_phenol = False
                has_alcohol = False

                for r_mol in reactant_mols:
                    if r_mol.HasSubstructMatch(phenol_pattern):
                        has_phenol = True
                    if r_mol.HasSubstructMatch(alcohol_pattern):
                        has_alcohol = True

                has_ether = product_mol.HasSubstructMatch(ether_pattern)

                # If we have a phenol and an alcohol as reactants, and an ether as product,
                # it's likely an O-alkylation
                if has_phenol and has_alcohol and has_ether:
                    has_o_alkylation = True
                    print(f"O-alkylation detected: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_o_alkylation

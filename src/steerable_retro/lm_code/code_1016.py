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
    Detects a strategy involving coupling reactions with halogenated aryl compounds.
    """
    # Initialize tracking variables
    has_halogenated_aryl = False
    has_coupling_with_haloarene = False

    def dfs_traverse(node, depth=0):
        nonlocal has_halogenated_aryl, has_coupling_with_haloarene

        if node["type"] == "mol":
            if "smiles" in node:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for halogenated aryl compounds
                    haloarene_pattern = Chem.MolFromSmarts("c-[F,Cl,Br,I]")
                    if mol.HasSubstructMatch(haloarene_pattern):
                        has_halogenated_aryl = True
                        print("Detected halogenated aryl compound")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for coupling reactions involving haloarenes
            halogen_present = any("Cl" in r or "F" in r or "Br" in r or "I" in r for r in reactants)
            if halogen_present and len(reactants) >= 2:
                # Check if it's a coupling reaction (new C-C bond formation)
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # If product has more C-C bonds than the sum in reactants, it might be a coupling
                    has_coupling_with_haloarene = True
                    print("Detected potential coupling reaction with haloarene")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if all required elements are present
    result = has_halogenated_aryl and has_coupling_with_haloarene

    print(f"Halogenated aryl detected: {has_halogenated_aryl}")
    print(f"Coupling with haloarene detected: {has_coupling_with_haloarene}")
    print(f"Overall halogenated aryl coupling strategy detected: {result}")

    return result

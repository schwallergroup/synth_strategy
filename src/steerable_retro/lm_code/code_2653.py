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
    This function detects if the route contains a TBDMS protection and deprotection sequence.
    TBDMS (tert-butyldimethylsilyl) is a common protecting group for alcohols.
    """
    has_protected = False
    has_deprotection = False

    def dfs_traverse(node):
        nonlocal has_protected, has_deprotection

        if node["type"] == "mol" and node.get("smiles"):
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for TBDMS protected alcohol
                tbdms_pattern = Chem.MolFromSmarts("[#14]([#6])([#6])([#6])[#8]-[#6]")
                if mol.HasSubstructMatch(tbdms_pattern):
                    has_protected = True
                    print(f"Found TBDMS protected intermediate: {node['smiles']}")

        elif node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for TBDMS deprotection (protected alcohol in reactants, free alcohol in product)
            reactant_has_tbdms = False
            for reactant in reactants:
                if reactant:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[#14]([#6])([#6])([#6])[#8]-[#6]")
                    ):
                        reactant_has_tbdms = True

            product_mol = Chem.MolFromSmiles(product)
            product_has_oh = product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts("[#8;H1]-[#6]")
            )

            if reactant_has_tbdms and product_has_oh:
                has_deprotection = True
                print(f"Found TBDMS deprotection: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Protection-deprotection sequence: Protected={has_protected}, Deprotected={has_deprotection}"
    )
    return has_protected and has_deprotection

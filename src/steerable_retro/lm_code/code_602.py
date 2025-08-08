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
    This function detects if the synthetic route involves methoxy deprotection
    to reveal a phenol group.
    """
    # Track if we found methoxy deprotection
    deprotection_found = False

    def dfs_traverse(node):
        nonlocal deprotection_found

        # Only process reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert SMILES to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and all(r is not None for r in reactants):
                # Check if any reactant contains methoxy on aromatic ring
                methoxy_pattern = Chem.MolFromSmarts("c-[O]-[CH3]")

                # Check if product contains phenol
                phenol_pattern = Chem.MolFromSmarts("c[OH]")

                has_methoxy = any(r.HasSubstructMatch(methoxy_pattern) for r in reactants)
                has_phenol = product.HasSubstructMatch(phenol_pattern) if product else False

                # If this is a methoxy deprotection reaction
                if has_methoxy and has_phenol:
                    print("Found methoxy deprotection to phenol")
                    deprotection_found = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return deprotection_found

#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects a synthetic strategy involving O-alkylation of a phenol.
    """
    phenol_alkylation_detected = False

    def dfs_traverse(node):
        nonlocal phenol_alkylation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                try:
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                    # Check for phenol alkylation pattern
                    phenol_pattern = Chem.MolFromSmarts("[c][OH]")
                    alkylated_phenol_pattern = Chem.MolFromSmarts("[c][O][C]")

                    if product_mol and phenol_pattern and alkylated_phenol_pattern:
                        if product_mol.HasSubstructMatch(alkylated_phenol_pattern):
                            # Check if any reactant has phenol
                            for r_mol in reactant_mols:
                                if r_mol and r_mol.HasSubstructMatch(phenol_pattern):
                                    phenol_alkylation_detected = True
                                    print("Phenol alkylation detected")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return phenol_alkylation_detected

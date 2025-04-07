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
    This function detects a synthetic strategy involving nucleophilic aromatic substitution
    to introduce an amine group, replacing a halogen.
    """
    snar_detected = False

    def dfs_traverse(node):
        nonlocal snar_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if all(reactants_mols) and product_mol:
                    # Check for halogen on aromatic in reactants
                    halo_aromatic_pattern = Chem.MolFromSmarts("[c][F,Cl,Br,I]")

                    # Check for amine nucleophile in reactants
                    amine_pattern = Chem.MolFromSmarts("[NH2][C,c]")

                    # Check for C-N bond in product where halogen was
                    cn_aromatic_pattern = Chem.MolFromSmarts("[c][N]")

                    reactants_with_halo = any(
                        [mol.HasSubstructMatch(halo_aromatic_pattern) for mol in reactants_mols]
                    )
                    reactants_with_amine = any(
                        [mol.HasSubstructMatch(amine_pattern) for mol in reactants_mols]
                    )
                    product_with_cn = product_mol.HasSubstructMatch(cn_aromatic_pattern)

                    if reactants_with_halo and reactants_with_amine and product_with_cn:
                        print("SNAr amine introduction detected")
                        snar_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return snar_detected

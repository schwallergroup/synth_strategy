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
    This function detects if the synthesis involves an SNAr reaction with a
    halogenated aromatic compound.
    """
    snar_detected = False

    def dfs_traverse(node):
        nonlocal snar_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Check for SNAr pattern: halogenated aromatic + amine → C-N bond formation
                    halo_aromatic_pattern = Chem.MolFromSmarts("c-[F,Cl,Br,I]")
                    amine_pattern = Chem.MolFromSmarts("[#7;H2]-c")
                    c_n_bond_pattern = Chem.MolFromSmarts("c-[#7]-c")

                    reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if reactants_mol and product_mol:
                        if (
                            reactants_mol.HasSubstructMatch(halo_aromatic_pattern)
                            and reactants_mol.HasSubstructMatch(amine_pattern)
                            and product_mol.HasSubstructMatch(c_n_bond_pattern)
                        ):
                            print("Detected SNAr with halogenated aromatic")
                            snar_detected = True
                except Exception as e:
                    print(f"Error in SNAr analysis: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return snar_detected

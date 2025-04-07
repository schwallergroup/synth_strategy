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
    This function detects if the synthesis involves a nucleophilic substitution
    between a thiol and a benzyl chloride.
    """
    thiol_benzyl_cl_coupling = False

    def dfs_traverse(node):
        nonlocal thiol_benzyl_cl_coupling

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for thiol pattern in reactants
            thiol_pattern = Chem.MolFromSmarts("[SH]")
            benzyl_cl_pattern = Chem.MolFromSmarts("[#6][CH2]Cl")

            thiol_present = False
            benzyl_cl_present = False

            for r_smi in reactants_smiles:
                try:
                    r_mol = Chem.MolFromSmiles(r_smi)
                    if r_mol:
                        if r_mol.HasSubstructMatch(thiol_pattern):
                            thiol_present = True
                        if r_mol.HasSubstructMatch(benzyl_cl_pattern):
                            benzyl_cl_present = True
                except:
                    continue

            # Check if product has C-S bond where the thiol and benzyl chloride connected
            if thiol_present and benzyl_cl_present:
                try:
                    p_mol = Chem.MolFromSmiles(product_smiles)
                    cs_bond_pattern = Chem.MolFromSmarts("[#6][CH2][S][#6]")
                    if p_mol and p_mol.HasSubstructMatch(cs_bond_pattern):
                        thiol_benzyl_cl_coupling = True
                        print("Detected thiol-benzyl chloride coupling")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return thiol_benzyl_cl_coupling

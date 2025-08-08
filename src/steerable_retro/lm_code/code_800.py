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
    This function detects a C-C bond formation via coupling reaction,
    particularly focusing on alkyne-aryl coupling.
    """
    cc_coupling_detected = False

    def dfs_traverse(node):
        nonlocal cc_coupling_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for terminal alkyne in reactants
            terminal_alkyne_present = False
            aryl_halide_present = False

            for r_smi in reactants_smiles:
                try:
                    r_mol = Chem.MolFromSmiles(r_smi)
                    if r_mol:
                        if r_mol.HasSubstructMatch(Chem.MolFromSmarts("C#[CH]")):
                            terminal_alkyne_present = True
                        if r_mol.HasSubstructMatch(Chem.MolFromSmarts("c[Br,I,Cl]")):
                            aryl_halide_present = True
                except:
                    continue

            # Check if product has aryl-alkyne bond
            if terminal_alkyne_present and aryl_halide_present:
                try:
                    p_mol = Chem.MolFromSmiles(product_smiles)
                    if p_mol and p_mol.HasSubstructMatch(Chem.MolFromSmarts("c[C]#[C]")):
                        cc_coupling_detected = True
                        print("C-C coupling (aryl-alkyne) detected")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return cc_coupling_detected

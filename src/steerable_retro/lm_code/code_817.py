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
    This function detects SNAr (nucleophilic aromatic substitution) reactions in the synthetic route.
    """
    snar_detected = False

    def dfs_traverse(node):
        nonlocal snar_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for aromatic halide and nucleophile in reactants, and substituted product
            aromatic_halide_pattern = Chem.MolFromSmarts("c[Cl,Br,F,I]")
            phenol_pattern = Chem.MolFromSmarts("c[OH]")
            ether_pattern = Chem.MolFromSmarts("cOc")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(r for r in reactant_mols):
                has_aromatic_halide = any(
                    r and r.HasSubstructMatch(aromatic_halide_pattern) for r in reactant_mols
                )
                has_phenol = any(r and r.HasSubstructMatch(phenol_pattern) for r in reactant_mols)
                has_ether = product_mol.HasSubstructMatch(ether_pattern)

                if has_aromatic_halide and has_phenol and has_ether:
                    print("Detected SNAr reaction")
                    snar_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return snar_detected

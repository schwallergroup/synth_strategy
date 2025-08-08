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
    This function detects a synthetic strategy involving C-N bond formation between
    heterocyclic compounds, particularly focusing on aryl-amine couplings.
    """
    # Track if we found the pattern
    found_cn_coupling = False

    def dfs_traverse(node):
        nonlocal found_cn_coupling

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for aryl-amine coupling
            # Look for aryl halide pattern in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,Cl,I,F]")
            # Look for amine pattern in reactants
            amine_pattern = Chem.MolFromSmarts("[NH2]")
            # Look for aryl-amine pattern in product
            aryl_amine_pattern = Chem.MolFromSmarts("[c][NH][c]")

            reactant_has_aryl_halide = False
            reactant_has_amine = False

            for r_smiles in reactants_smiles:
                try:
                    r_mol = Chem.MolFromSmiles(r_smiles)
                    if r_mol:
                        if r_mol.HasSubstructMatch(aryl_halide_pattern):
                            reactant_has_aryl_halide = True
                        if r_mol.HasSubstructMatch(amine_pattern):
                            reactant_has_amine = True
                except:
                    continue

            product_has_aryl_amine = False
            try:
                p_mol = Chem.MolFromSmiles(product_smiles)
                if p_mol and p_mol.HasSubstructMatch(aryl_amine_pattern):
                    product_has_aryl_amine = True
            except:
                pass

            # Check for heterocycles in reactants
            pyrazole_pattern = Chem.MolFromSmarts("[#6]1[#6][#7][#7][#6]1")
            pyrimidine_pattern = Chem.MolFromSmarts("[#6]1[#6][#7][#6][#7][#6]1")

            reactant_has_heterocycle = False
            for r_smiles in reactants_smiles:
                try:
                    r_mol = Chem.MolFromSmiles(r_smiles)
                    if r_mol:
                        if r_mol.HasSubstructMatch(pyrazole_pattern) or r_mol.HasSubstructMatch(
                            pyrimidine_pattern
                        ):
                            reactant_has_heterocycle = True
                            break
                except:
                    continue

            if (
                reactant_has_aryl_halide
                and reactant_has_amine
                and product_has_aryl_amine
                and reactant_has_heterocycle
            ):
                found_cn_coupling = True
                print("Found heterocycle C-N coupling")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_cn_coupling

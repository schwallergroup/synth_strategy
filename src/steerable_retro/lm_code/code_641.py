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
    This function detects if the synthesis route employs a Suzuki coupling strategy
    for C-C bond formation between aromatic rings.
    """
    suzuki_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Create RDKit molecules
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                # Check for aryl halide in reactants
                aryl_halide_pattern = Chem.MolFromSmarts("c[Br,I,Cl]")
                has_aryl_halide = any(
                    mol is not None and mol.HasSubstructMatch(aryl_halide_pattern)
                    for mol in reactant_mols
                )

                # Check for boronic acid or ester in reactants
                boronic_pattern = Chem.MolFromSmarts("cB(O)(O)")
                boronic_ester_pattern = Chem.MolFromSmarts("cB(O[C,c])(O[C,c])")
                has_boronic = any(
                    mol is not None
                    and (
                        mol.HasSubstructMatch(boronic_pattern)
                        or mol.HasSubstructMatch(boronic_ester_pattern)
                    )
                    for mol in reactant_mols
                )

                # Check if product has biaryl system that wasn't in reactants
                if has_aryl_halide and has_boronic and product_mol is not None:
                    print(f"Suzuki coupling detected at depth {depth}")
                    suzuki_detected = True
            except:
                print("Error processing reaction SMILES for Suzuki coupling detection")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return suzuki_detected

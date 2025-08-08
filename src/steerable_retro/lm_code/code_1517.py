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
    This function detects the use of Suzuki coupling (aryl halide + boronic acid)
    to form C-C bonds between aromatic rings.
    """
    # Track if we found a Suzuki coupling
    found_suzuki_coupling = False

    # Define SMARTS patterns
    aryl_halide_pattern = Chem.MolFromSmarts("c[Br,I,Cl]")
    boronic_acid_pattern = Chem.MolFromSmarts("cB(O)O")
    biaryl_pattern = Chem.MolFromSmarts("c:c")

    def dfs_traverse(node):
        nonlocal found_suzuki_coupling

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Convert to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                # Check if reactants include aryl halide and boronic acid
                has_aryl_halide = any(
                    mol is not None and mol.HasSubstructMatch(aryl_halide_pattern)
                    for mol in reactant_mols
                )
                has_boronic_acid = any(
                    mol is not None and mol.HasSubstructMatch(boronic_acid_pattern)
                    for mol in reactant_mols
                )

                # Check if product has a biaryl bond
                has_biaryl = product_mol is not None and product_mol.HasSubstructMatch(
                    biaryl_pattern
                )

                # If we have both required reactants and the biaryl product
                if has_aryl_halide and has_boronic_acid and has_biaryl:
                    found_suzuki_coupling = True
                    print("Found Suzuki coupling reaction")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return found_suzuki_coupling

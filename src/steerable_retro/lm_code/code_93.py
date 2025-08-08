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
    This function detects SNAr reactions for C-N bond formation,
    typically involving nucleophilic substitution on aromatic rings with electron-withdrawing groups.
    """
    has_snar_reaction = False

    def dfs_traverse(node, depth=0):
        nonlocal has_snar_reaction

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            if product_smiles and Chem.MolFromSmiles(product_smiles):
                product_mol = Chem.MolFromSmiles(product_smiles)

                # Check for amine reactant
                has_amine = False
                has_aryl_halide = False
                amine_reactant = None

                for reactant in reactants_smiles:
                    if reactant and Chem.MolFromSmiles(reactant):
                        reactant_mol = Chem.MolFromSmiles(reactant)

                        # Check for primary or secondary amine
                        amine_pattern = Chem.MolFromSmarts("[NH1,NH2]")
                        if reactant_mol.HasSubstructMatch(amine_pattern):
                            has_amine = True
                            amine_reactant = reactant
                            print(f"Found amine reactant in reaction at depth {depth}")

                        # Check for aryl halide
                        aryl_halide_pattern = Chem.MolFromSmarts("c[Cl,F,Br,I]")
                        if reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True
                            print(f"Found aryl halide in reaction at depth {depth}")

                # Check for C-N bond formation
                if has_amine and has_aryl_halide and amine_reactant:
                    # Check if product has new C-N bond
                    # This is a simplification - would need more sophisticated atom mapping
                    c_n_bond_pattern = Chem.MolFromSmarts("c[NH0,NH1]")
                    if product_mol.HasSubstructMatch(c_n_bond_pattern):
                        has_snar_reaction = True
                        print(f"Detected SNAr C-N bond formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_snar_reaction

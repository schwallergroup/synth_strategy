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
    This function detects a synthetic strategy involving multiple sequential C-N bond formations,
    including at least one N-alkylation and one reductive amination.
    """
    c_n_bond_formations = 0
    has_reductive_amination = False
    has_n_alkylation = False

    def dfs_traverse(node):
        nonlocal c_n_bond_formations, has_reductive_amination, has_n_alkylation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for C-N bond formation
            if product and all(r for r in reactants):
                # Look for aldehyde to amine conversion (reductive amination)
                aldehyde_pattern = Chem.MolFromSmarts("[C;H1]=O")
                amine_pattern = Chem.MolFromSmarts("[#6][N;!$(NC=O)]")

                has_aldehyde = any(r.HasSubstructMatch(aldehyde_pattern) for r in reactants if r)
                has_amine_product = product.HasSubstructMatch(amine_pattern)

                if has_aldehyde and has_amine_product:
                    has_reductive_amination = True
                    c_n_bond_formations += 1
                    print("Detected reductive amination")

                # Check for N-alkylation (alkyl halide + amine)
                alkyl_halide_pattern = Chem.MolFromSmarts("[#6][C;!$(C=O)][F,Cl,Br,I]")
                amine_reactant_pattern = Chem.MolFromSmarts("[N;!$(NC=O)]")

                has_alkyl_halide = any(
                    r.HasSubstructMatch(alkyl_halide_pattern) for r in reactants if r
                )
                has_amine_reactant = any(
                    r.HasSubstructMatch(amine_reactant_pattern) for r in reactants if r
                )

                if has_alkyl_halide and has_amine_reactant and has_amine_product:
                    has_n_alkylation = True
                    c_n_bond_formations += 1
                    print("Detected N-alkylation")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we have multiple C-N bond formations including both types
    return c_n_bond_formations >= 2 and has_reductive_amination and has_n_alkylation

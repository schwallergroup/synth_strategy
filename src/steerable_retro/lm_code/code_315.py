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
    This function detects if the synthetic route employs a strategy involving
    protection and deprotection of phenol groups.
    """
    # Initialize tracking variables
    phenol_protection_steps = []
    phenol_deprotection_steps = []

    # Define SMARTS patterns
    phenol_pattern = Chem.MolFromSmarts("c[OH]")
    protected_phenol_pattern = Chem.MolFromSmarts("cO[!H]")  # Phenol with non-hydrogen substituent

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for phenol protection (phenol in reactants, protected phenol in product)
                has_phenol_reactant = False
                for reactant_smiles in reactants_smiles:
                    try:
                        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                        if reactant_mol and reactant_mol.HasSubstructMatch(phenol_pattern):
                            has_phenol_reactant = True
                            break
                    except:
                        print(f"Error processing reactant SMILES: {reactant_smiles}")

                has_protected_phenol_product = False
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(protected_phenol_pattern):
                        has_protected_phenol_product = True
                except:
                    print(f"Error processing product SMILES: {product_smiles}")

                if has_phenol_reactant and has_protected_phenol_product:
                    print(f"Detected phenol protection at depth {depth}")
                    phenol_protection_steps.append(depth)

                # Check for phenol deprotection (protected phenol in reactants, phenol in product)
                has_protected_phenol_reactant = False
                for reactant_smiles in reactants_smiles:
                    try:
                        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            protected_phenol_pattern
                        ):
                            has_protected_phenol_reactant = True
                            break
                    except:
                        print(f"Error processing reactant SMILES: {reactant_smiles}")

                has_phenol_product = False
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(phenol_pattern):
                        has_phenol_product = True
                except:
                    print(f"Error processing product SMILES: {product_smiles}")

                if has_protected_phenol_reactant and has_phenol_product:
                    print(f"Detected phenol deprotection at depth {depth}")
                    phenol_deprotection_steps.append(depth)

        # Continue traversing children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if both protection and deprotection occurred
    has_protection_deprotection = (
        len(phenol_protection_steps) > 0 and len(phenol_deprotection_steps) > 0
    )

    if has_protection_deprotection:
        print("Detected phenol protection-deprotection strategy")

    return has_protection_deprotection

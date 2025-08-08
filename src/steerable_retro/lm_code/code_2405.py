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
    Detects a synthetic strategy involving exchange between different halogenated functional groups.
    """
    has_halogen_exchange = False

    # SMARTS patterns
    chloride_pattern = Chem.MolFromSmarts("[#6]-[Cl]")
    bromide_pattern = Chem.MolFromSmarts("[#6]-[Br]")
    alcohol_pattern = Chem.MolFromSmarts("[#6]-[OH]")

    def dfs_traverse(node):
        nonlocal has_halogen_exchange

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product_mol and reactant_mols:
                # Check for various halogen exchanges
                reactant_has_chloride = any(
                    r and r.HasSubstructMatch(chloride_pattern) for r in reactant_mols
                )
                reactant_has_bromide = any(
                    r and r.HasSubstructMatch(bromide_pattern) for r in reactant_mols
                )
                reactant_has_alcohol = any(
                    r and r.HasSubstructMatch(alcohol_pattern) for r in reactant_mols
                )

                product_has_chloride = product_mol and product_mol.HasSubstructMatch(
                    chloride_pattern
                )
                product_has_bromide = product_mol and product_mol.HasSubstructMatch(bromide_pattern)
                product_has_alcohol = product_mol and product_mol.HasSubstructMatch(alcohol_pattern)

                # Check for specific exchanges
                if (
                    (reactant_has_chloride and (product_has_alcohol or product_has_bromide))
                    or (reactant_has_bromide and (product_has_alcohol or product_has_chloride))
                    or (reactant_has_alcohol and (product_has_chloride or product_has_bromide))
                ):
                    has_halogen_exchange = True
                    print("Detected halogen exchange")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Halogen exchange strategy detected: {has_halogen_exchange}")
    return has_halogen_exchange

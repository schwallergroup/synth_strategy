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
    This function detects a synthetic strategy involving heterocyclic ring closure
    via C-N bond formation to construct a nitrogen-containing ring system.
    """
    # Initialize tracking variables
    has_cn_bond_formation = False
    has_heterocycle_formation = False

    def dfs_traverse(node):
        nonlocal has_cn_bond_formation, has_heterocycle_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product_mol and len(reactant_mols) > 1:
                # Check for C-N bond formation
                # This is a simplified approach - in practice, you'd need atom mapping
                # to accurately track bond formation

                # Check if product has more rings than reactants
                product_ring_count = sum(1 for ring in Chem.GetSSSR(product_mol))
                reactants_ring_count = sum(
                    sum(1 for ring in Chem.GetSSSR(r)) for r in reactant_mols
                )

                if product_ring_count > reactants_ring_count:
                    # Check if the product contains a nitrogen heterocycle
                    n_heterocycle_pattern = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#7]~[#6]1")
                    if product_mol.HasSubstructMatch(n_heterocycle_pattern):
                        has_heterocycle_formation = True
                        has_cn_bond_formation = True
                        print("Detected heterocyclic ring formation via C-N bond")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = has_cn_bond_formation and has_heterocycle_formation

    if strategy_present:
        print("Detected heterocyclic ring closure via C-N bond formation")

    return strategy_present

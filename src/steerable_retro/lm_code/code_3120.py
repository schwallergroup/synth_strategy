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
    Detects a strategy where an amino group is introduced in the late stage of synthesis
    (at depth 0 or 1).
    """
    late_stage_amination_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amination_detected

        if node["type"] == "reaction" and depth <= 1:  # Only check reactions at depth 0 or 1
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]

                # Check for amino group in product
                amino_pattern = Chem.MolFromSmarts("[NH2]")
                if product_mol and product_mol.HasSubstructMatch(amino_pattern):
                    # Check if amino group is not present in reactants
                    amino_in_reactants = False
                    for r_mol in reactant_mols:
                        if r_mol and r_mol.HasSubstructMatch(amino_pattern):
                            amino_in_reactants = True
                            break

                    if not amino_in_reactants:
                        print(f"Detected late-stage amination at depth {depth}")
                        late_stage_amination_detected = True

            except Exception as e:
                print(f"Error in SMILES processing: {e}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return late_stage_amination_detected

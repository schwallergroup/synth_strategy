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
    This function detects a synthetic strategy involving the formation of a ketone
    that bridges two aromatic systems.
    """
    # Track if we find a ketone connecting two aromatic rings
    has_aromatic_carbonyl_bridge = False

    # Define SMARTS patterns
    diaryl_ketone_pattern = Chem.MolFromSmarts("c-C(=O)-c")

    def dfs_traverse(node):
        nonlocal has_aromatic_carbonyl_bridge

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Check if product contains a ketone connecting two aromatic rings
                    if product_mol and product_mol.HasSubstructMatch(diaryl_ketone_pattern):
                        # Check if this is a new formation (not present in reactants)
                        reactant_has_pattern = False
                        for r_smiles in reactants_smiles:
                            if r_smiles:
                                r_mol = Chem.MolFromSmiles(r_smiles)
                                if r_mol and r_mol.HasSubstructMatch(diaryl_ketone_pattern):
                                    reactant_has_pattern = True
                                    break

                        if not reactant_has_pattern:
                            print("Detected formation of aromatic-carbonyl-aromatic bridge")
                            has_aromatic_carbonyl_bridge = True

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Aromatic carbonyl bridge strategy detected: {has_aromatic_carbonyl_bridge}")
    return has_aromatic_carbonyl_bridge

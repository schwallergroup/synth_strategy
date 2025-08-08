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
    This function detects a synthetic strategy involving heterocyclic ring opening
    in the early stages of synthesis.
    """
    ring_opening_detected = False

    def dfs_traverse(node):
        nonlocal ring_opening_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert SMILES to molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if not all(reactant_mols) or not product_mol:
                print("Warning: Could not parse all molecules in reaction")
                return

            # Count rings in reactants and product
            reactant_ring_counts = [rdMolDescriptors.CalcNumRings(mol) for mol in reactant_mols]
            product_ring_count = rdMolDescriptors.CalcNumRings(product_mol)

            total_reactant_rings = sum(reactant_ring_counts)

            # Check if product has fewer rings than reactants (ring opening)
            if product_ring_count < total_reactant_rings:
                # Check if the molecules contain nitrogen (heterocycles)
                nitrogen_in_reactants = any(
                    mol.GetSubstructMatches(Chem.MolFromSmarts("[#7]")) for mol in reactant_mols
                )
                nitrogen_in_product = bool(
                    product_mol.GetSubstructMatches(Chem.MolFromSmarts("[#7]"))
                )

                if nitrogen_in_reactants and nitrogen_in_product:
                    print(f"Heterocycle ring opening detected: {rsmi}")
                    ring_opening_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return ring_opening_detected

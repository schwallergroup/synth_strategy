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
    This function detects a synthetic strategy involving the formation of a cyclopentenone ring
    from an acyclic precursor.
    """
    cyclopentenone_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal cyclopentenone_formation_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert SMILES to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product_mol and reactant_mols:
                # Check for cyclopentenone formation
                cyclopentenone_pattern = Chem.MolFromSmarts("O=C1CCC=C1")

                if product_mol.HasSubstructMatch(cyclopentenone_pattern) and not any(
                    r.HasSubstructMatch(cyclopentenone_pattern) for r in reactant_mols
                ):
                    print(f"Cyclopentenone formation detected at depth {depth}")
                    cyclopentenone_formation_found = True

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Cyclopentenone formation strategy detected: {cyclopentenone_formation_found}")
    return cyclopentenone_formation_found

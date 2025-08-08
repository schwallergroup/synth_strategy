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
    This function detects hydrazide formation from ester precursors.
    """
    hydrazide_formation_found = False

    def dfs_traverse(node):
        nonlocal hydrazide_formation_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                products_smiles = rsmi.split(">")[-1]

                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".") if r]
                products = [Chem.MolFromSmiles(p) for p in products_smiles.split(".") if p]

                # Ester pattern
                ester_patt = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])")

                # Hydrazide pattern
                hydrazide_patt = Chem.MolFromSmarts("[#6](=[#8])-[#7]-[#7]")

                # Hydrazine pattern
                hydrazine_patt = Chem.MolFromSmarts("[#7]-[#7]")

                # Check for ester in reactants and hydrazide in products
                # Also check for hydrazine in reactants
                if (
                    any(mol.HasSubstructMatch(ester_patt) for mol in reactants if mol)
                    and any(mol.HasSubstructMatch(hydrazine_patt) for mol in reactants if mol)
                    and any(mol.HasSubstructMatch(hydrazide_patt) for mol in products if mol)
                ):
                    hydrazide_formation_found = True
                    print("Found hydrazide formation from ester")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return hydrazide_formation_found

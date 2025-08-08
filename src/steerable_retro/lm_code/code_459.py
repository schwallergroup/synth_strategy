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
    Detects a convergent synthesis strategy where a cyano-containing fragment
    is combined with another fragment through amide bond formation.
    """
    found_convergent_step = False
    found_cyano_fragment = False

    def dfs_traverse(node, depth=0):
        nonlocal found_convergent_step, found_cyano_fragment

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is a convergent step (multiple reactants)
            if len(reactants_smiles) > 1:
                # Convert to RDKit molecules
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]

                # Check for cyano group in reactants
                cyano_pattern = Chem.MolFromSmarts("[C]#[N]")
                has_cyano = any(r and r.HasSubstructMatch(cyano_pattern) for r in reactants if r)

                # Check for amide formation
                product = Chem.MolFromSmiles(product_smiles) if product_smiles else None
                amide_pattern = Chem.MolFromSmarts("[C](=[O])[NH][c]")
                has_amide = product and product.HasSubstructMatch(amide_pattern)

                if has_cyano and has_amide:
                    found_convergent_step = True
                    found_cyano_fragment = True
                    print("Found convergent synthesis with cyano fragment")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_convergent_step and found_cyano_fragment

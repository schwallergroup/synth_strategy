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
    This function detects if the synthetic route involves olefination chemistry
    to form a C=C bond adjacent to a nitrile group.
    """
    nitrile_olefination_found = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_olefination_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                products = rsmi.split(">")[-1]

                # Check for ketone/aldehyde in reactants
                reactant_mol = Chem.MolFromSmiles(reactants)
                carbonyl_pattern = Chem.MolFromSmarts("[#6][#6](=[O])[#6]")

                # Check for α,β-unsaturated nitrile in products
                product_mol = Chem.MolFromSmiles(products)
                unsaturated_nitrile_pattern = Chem.MolFromSmarts("[#6]=[#6][#6]#[N]")

                if (
                    reactant_mol
                    and carbonyl_pattern
                    and reactant_mol.HasSubstructMatch(carbonyl_pattern)
                    and product_mol
                    and unsaturated_nitrile_pattern
                    and product_mol.HasSubstructMatch(unsaturated_nitrile_pattern)
                ):
                    print(f"Found nitrile olefination strategy at depth {depth}")
                    nitrile_olefination_found = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return nitrile_olefination_found

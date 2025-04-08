#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects a late-stage reductive amination strategy.
    It looks for a reaction at depth 0 or 1 that combines a carbonyl compound with an amine.
    """
    reductive_amination_found = False

    def dfs_traverse(node, depth=0):
        nonlocal reductive_amination_found

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Check if reactants contain both carbonyl and amine groups
                reactants = reactants_str.split(".")

                carbonyl_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#6]")
                amine_pattern = Chem.MolFromSmarts("[#7H2]-[#6]")

                has_carbonyl = False
                has_amine = False

                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(carbonyl_pattern):
                            has_carbonyl = True
                        if mol and mol.HasSubstructMatch(amine_pattern):
                            has_amine = True
                    except:
                        continue

                # Check if product has C-N bond that wasn't in reactants
                if has_carbonyl and has_amine:
                    try:
                        product_mol = Chem.MolFromSmiles(product_str)
                        cn_bond_pattern = Chem.MolFromSmarts("[#6]-[#7]-[#6]")
                        if product_mol and product_mol.HasSubstructMatch(cn_bond_pattern):
                            print("Late-stage reductive amination detected at depth", depth)
                            reductive_amination_found = True
                    except:
                        pass

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return reductive_amination_found

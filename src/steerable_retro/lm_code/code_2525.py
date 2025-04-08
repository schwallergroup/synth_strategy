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
    Detects if the synthesis uses multiple N-acylation reactions.
    """
    acylation_count = 0

    def dfs_traverse(node):
        nonlocal acylation_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an acylation reaction
                acyl_pattern = Chem.MolFromSmarts("[#7:1][C:2](=[O:3])[#6:4]")
                prod_mol = Chem.MolFromSmiles(product)

                if prod_mol and prod_mol.HasSubstructMatch(acyl_pattern):
                    # Check if reactants include an acid chloride or similar acylating agent
                    acyl_agent_pattern = Chem.MolFromSmarts("[Cl,Br,I,OH][C](=[O])[#6]")
                    amine_pattern = Chem.MolFromSmarts("[#7;!$(N~[!#6])]")

                    has_acyl_agent = False
                    has_amine = False

                    for reactant in reactants:
                        r_mol = Chem.MolFromSmiles(reactant)
                        if r_mol:
                            if r_mol.HasSubstructMatch(acyl_agent_pattern):
                                has_acyl_agent = True
                            if r_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True

                    if has_acyl_agent and has_amine:
                        acylation_count += 1
                        print(f"Found acylation reaction #{acylation_count}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return acylation_count >= 2

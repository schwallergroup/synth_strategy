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
    Detects if the synthesis route involves N-alkylation with an alpha-bromoester
    or similar activated alkylating agent.
    """

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Look for alpha-bromoester pattern in reactants
                    alpha_bromoester_pattern = Chem.MolFromSmarts(
                        "[Br,Cl,I][#6][#6](=[#8])[#8][#6]"
                    )
                    amine_pattern = Chem.MolFromSmarts(
                        "[#7;!$(N~[!#6]);!$(N~[#6]=[#8])]"
                    )  # Amine not amide

                    # Check if any reactant contains alpha-bromoester
                    has_alpha_bromoester = False
                    has_amine = False

                    for r in reactants:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol is None:
                            continue

                        if r_mol.HasSubstructMatch(alpha_bromoester_pattern):
                            has_alpha_bromoester = True
                            print(f"Alpha-bromoester found in reactant: {r}")

                        if r_mol.HasSubstructMatch(amine_pattern):
                            has_amine = True
                            print(f"Amine found in reactant: {r}")

                    # Check if product has a new C-N bond
                    if has_alpha_bromoester and has_amine:
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol is not None:
                            # This is a simplified check - in a real implementation,
                            # you would need to track atom mappings to confirm the specific bond formation
                            print(f"N-alkylation with alpha-bromoester detected at depth {depth}")
                            return True
                except:
                    pass

        for child in node.get("children", []):
            result = dfs_traverse(child, depth + 1)
            if result:
                return True

        return False

    result = dfs_traverse(route)
    if not result:
        print("No N-alkylation with alpha-bromoester detected")
    return result

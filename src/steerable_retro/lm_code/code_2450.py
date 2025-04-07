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
    This function detects if the synthesis involves an O-demethylation step
    (conversion of methoxy to hydroxyl).
    """
    has_demethylation = False

    def dfs_traverse(node):
        nonlocal has_demethylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Create molecules
                reactants_mol = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
                product_mol = Chem.MolFromSmiles(product_part)

                if product_mol:
                    # Define patterns
                    methoxy_pattern = Chem.MolFromSmarts("[c][O][C;H3]")
                    phenol_pattern = Chem.MolFromSmarts("[c][O;H1]")

                    # Check if product has phenol but not methoxy
                    if product_mol.HasSubstructMatch(phenol_pattern):
                        # Check if any reactant has methoxy
                        for r_mol in reactants_mol:
                            if r_mol and r_mol.HasSubstructMatch(methoxy_pattern):
                                has_demethylation = True
                                print(f"Found O-demethylation step: {rsmi}")
                                break

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Synthesis {'includes' if has_demethylation else 'does not include'} O-demethylation")
    return has_demethylation

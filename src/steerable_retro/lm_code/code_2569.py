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
    Detects a strategy involving multiple heteroatom alkylations (O-alkylation and N-alkylation).
    """
    o_alkylation_count = 0
    n_alkylation_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal o_alkylation_count, n_alkylation_count

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                try:
                    # Check for O-alkylation
                    if "[OH]" in reactants_str and "[O][C]" in product_str:
                        o_alkylation_count += 1
                        print(f"Detected O-alkylation at depth {depth}")

                    # Check for N-alkylation
                    if "[NH]" in reactants_str and "[N][C]" in product_str:
                        n_alkylation_count += 1
                        print(f"Detected N-alkylation at depth {depth}")

                    # Alternative detection using SMARTS
                    reactants = reactants_str.split(".")
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        product_mol = Chem.MolFromSmiles(product_str)

                        if reactant_mol and product_mol:
                            # O-alkylation patterns
                            phenol_pattern = Chem.MolFromSmarts("[OH]c")
                            alkoxy_pattern = Chem.MolFromSmarts("[O][CH2]")

                            # N-alkylation patterns
                            amine_pattern = Chem.MolFromSmarts("[NH]")
                            alkyl_amine_pattern = Chem.MolFromSmarts("[N][CH2]")

                            if reactant_mol.HasSubstructMatch(
                                phenol_pattern
                            ) and product_mol.HasSubstructMatch(alkoxy_pattern):
                                o_alkylation_count += 1
                                print(f"Detected O-alkylation at depth {depth} using SMARTS")

                            if reactant_mol.HasSubstructMatch(
                                amine_pattern
                            ) and product_mol.HasSubstructMatch(alkyl_amine_pattern):
                                n_alkylation_count += 1
                                print(f"Detected N-alkylation at depth {depth} using SMARTS")
                except:
                    print(f"Error processing alkylation detection at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have multiple heteroatom alkylations
    total_alkylations = o_alkylation_count + n_alkylation_count
    has_multiple_heteroatom_alkylations = total_alkylations >= 2

    print(
        f"Multiple heteroatom alkylation strategy detected: {has_multiple_heteroatom_alkylations}"
    )
    print(f"- O-alkylations: {o_alkylation_count}")
    print(f"- N-alkylations: {n_alkylation_count}")
    print(f"- Total alkylations: {total_alkylations}")

    return has_multiple_heteroatom_alkylations

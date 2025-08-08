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
    This function detects a synthetic strategy with multiple amide bond formations
    throughout the synthesis.
    """
    amide_formation_count = 0

    def dfs_traverse(node):
        nonlocal amide_formation_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide formation
            prod_mol = Chem.MolFromSmiles(product)
            if prod_mol:
                amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
                if prod_mol.HasSubstructMatch(amide_pattern):
                    # Check if reactants contain amine and carboxylic acid/acyl chloride
                    has_amine = False
                    has_carbonyl = False

                    for reactant in reactants:
                        react_mol = Chem.MolFromSmiles(reactant)
                        if react_mol:
                            amine_pattern = Chem.MolFromSmarts("[N;!$(N=*);!$(NC=O)]")
                            carbonyl_pattern = Chem.MolFromSmarts("[C](=[O])[O,Cl]")

                            if react_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True
                            if react_mol.HasSubstructMatch(carbonyl_pattern):
                                has_carbonyl = True

                    if has_amine and has_carbonyl:
                        amide_formation_count += 1
                        print(
                            f"Detected amide formation at depth {node['metadata'].get('depth', 'unknown')}"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if there are at least 2 amide formations
    return amide_formation_count >= 2

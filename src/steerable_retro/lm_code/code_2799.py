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
    Detects the use of reductive amination (aldehyde/ketone + amine → amine)
    as a key step in the synthetic route.
    """
    found_reductive_amination = False

    def dfs_traverse(node):
        nonlocal found_reductive_amination

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                if product_mol and all(reactant_mols):
                    # Check for primary or secondary amine in product
                    amine_pattern = Chem.MolFromSmarts("[C][N;!$(N=*);!$(N#*)]")

                    # Check for aldehyde or ketone in reactants
                    carbonyl_pattern = Chem.MolFromSmarts("[C](=[O])[C,H]")

                    # Check for ammonia or amine in reactants
                    nitrogen_source_pattern = Chem.MolFromSmarts("[N;!$(N=*);!$(N#*)]")

                    if (
                        product_mol.HasSubstructMatch(amine_pattern)
                        and any(r.HasSubstructMatch(carbonyl_pattern) for r in reactant_mols)
                        and any(r.HasSubstructMatch(nitrogen_source_pattern) for r in reactant_mols)
                    ):

                        # Additional check: product should have C-N bond where C was previously part of C=O
                        print("Found potential reductive amination")
                        found_reductive_amination = True
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_reductive_amination

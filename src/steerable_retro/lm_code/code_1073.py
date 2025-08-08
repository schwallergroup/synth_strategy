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
    Detects if the synthesis route involves nucleophilic aromatic substitution.
    Specifically looks for replacement of Cl-Ar by O-Ar.
    """
    snar_found = False

    def dfs_traverse(node, depth=0):
        nonlocal snar_found

        if node["type"] == "reaction" and depth <= 1:  # Focus on late-stage reactions
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Patterns for SNAr
                cl_ar_pattern = Chem.MolFromSmarts("[#17]-c:n")  # Chloro-heteroaromatic
                o_ar_pattern = Chem.MolFromSmarts("[#8]-c:n")  # Oxygen-heteroaromatic

                product_mol = Chem.MolFromSmiles(product)

                # Check if any reactant has Cl-Ar that's replaced by O-Ar in the product
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol or not product_mol:
                        continue

                    if (
                        reactant_mol.HasSubstructMatch(cl_ar_pattern)
                        and product_mol.HasSubstructMatch(o_ar_pattern)
                        and not reactant_mol.HasSubstructMatch(o_ar_pattern)
                    ):
                        print(f"Found nucleophilic aromatic substitution at depth {depth}")
                        snar_found = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return snar_found

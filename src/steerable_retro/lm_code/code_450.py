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
    Detects if the route involves oxidation of a methyl group to a carbonyl group.
    """
    found_oxidation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_oxidation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # This is a simplified check - in a real implementation,
                # we would need to track atom mappings to confirm the transformation
                methyl_pattern = Chem.MolFromSmarts("[CH3]")
                carbonyl_pattern = Chem.MolFromSmarts("[#6](=[#8])")

                # Check if reactant has methyl group
                has_methyl = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(methyl_pattern):
                        has_methyl = True
                        break

                # Check if product has carbonyl group
                product_mol = Chem.MolFromSmiles(product)
                has_carbonyl = product_mol and product_mol.HasSubstructMatch(carbonyl_pattern)

                if has_methyl and has_carbonyl:
                    found_oxidation = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Methyl to carbonyl oxidation detection: {found_oxidation}")
    return found_oxidation

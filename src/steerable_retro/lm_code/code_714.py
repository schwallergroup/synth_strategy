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
    Detects if the synthetic route employs a convergent synthesis approach
    where two complex fragments are joined via biaryl formation.
    """
    # Track if we found a convergent biaryl synthesis
    found_convergent_biaryl = False

    def dfs_traverse(node, depth=0):
        nonlocal found_convergent_biaryl

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Need at least 2 reactants for convergent synthesis
                if len(reactants) >= 2:
                    # Check if both reactants are complex (have more than 15 atoms)
                    complex_reactants = 0
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.GetNumAtoms() > 15:
                                complex_reactants += 1
                        except:
                            continue

                    # Check if product has a biaryl bond
                    try:
                        prod_mol = Chem.MolFromSmiles(product)
                        if prod_mol:
                            biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")
                            if prod_mol.HasSubstructMatch(biaryl_pattern):
                                # If we have at least 2 complex reactants and a biaryl in product
                                if complex_reactants >= 2:
                                    found_convergent_biaryl = True
                                    print(f"Found convergent biaryl synthesis at depth {depth}")
                    except:
                        pass

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_convergent_biaryl

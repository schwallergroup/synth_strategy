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
    Detects a convergent synthesis where two fragments are joined via an ether linkage.
    Specifically looks for a reaction that forms a C-O bond between an aromatic ring and an alkyl chain.
    """
    convergent_synthesis_detected = False

    def dfs_traverse(node):
        nonlocal convergent_synthesis_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Check if this is a reaction that forms an ether linkage
                reactants = reactants_part.split(".")

                # If we have multiple reactants
                if len(reactants) >= 2:
                    product_mol = Chem.MolFromSmiles(product_part)

                    # Check for ether linkage in product
                    ether_pattern = Chem.MolFromSmarts("c-O-C")
                    if product_mol and product_mol.HasSubstructMatch(ether_pattern):
                        # Check if reactants are complex enough (not just simple reagents)
                        complex_reactants = 0
                        for reactant in reactants:
                            r_mol = Chem.MolFromSmiles(reactant)
                            if (
                                r_mol and r_mol.GetNumAtoms() > 5
                            ):  # Arbitrary threshold for "complex"
                                complex_reactants += 1

                        if complex_reactants >= 2:
                            print("Detected convergent synthesis with ether linkage formation")
                            convergent_synthesis_detected = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return convergent_synthesis_detected

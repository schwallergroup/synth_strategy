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
    This function detects nitration of an aromatic ring in the synthesis route.
    """
    # Track if we found nitration reaction
    found_nitration = False

    def dfs_traverse(node):
        nonlocal found_nitration

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitration pattern
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol:
                        nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")
                        if prod_mol.HasSubstructMatch(nitro_pattern):
                            # Check if nitro group was added in this step
                            nitro_found_in_reactants = False
                            for reactant in reactants:
                                if "[N+](=O)[O-]" in reactant:
                                    nitro_found_in_reactants = True
                                    break

                            if not nitro_found_in_reactants:
                                found_nitration = True
                                print(f"Found nitration reaction: {rsmi}")
                except:
                    pass

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Nitro functionalization strategy detected: {found_nitration}")
    return found_nitration

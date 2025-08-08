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
    This function detects preservation of stereocenters throughout the synthesis.
    """
    has_stereocenter = False
    preserves_stereocenter = True

    def dfs_traverse(node):
        nonlocal has_stereocenter, preserves_stereocenter

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    product_stereocenters = Chem.FindMolChiralCenters(
                        product_mol, includeUnassigned=False
                    )

                    if product_stereocenters:
                        has_stereocenter = True

                        # Check if any reactant also has the stereocenter
                        reactant_has_stereocenter = False
                        for r in reactants:
                            if r:
                                reactant_mol = Chem.MolFromSmiles(r)
                                if reactant_mol:
                                    reactant_stereocenters = Chem.FindMolChiralCenters(
                                        reactant_mol, includeUnassigned=False
                                    )
                                    if reactant_stereocenters:
                                        reactant_has_stereocenter = True
                                        break

                        if not reactant_has_stereocenter:
                            # If product has stereocenter but reactants don't, it's not preservation
                            preserves_stereocenter = False
                except:
                    print("Error processing reaction SMILES:", rsmi)

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is true if there are stereocenters and they are preserved
    result = has_stereocenter and preserves_stereocenter

    if result:
        print("Stereocenter-preserving strategy detected")

    return result

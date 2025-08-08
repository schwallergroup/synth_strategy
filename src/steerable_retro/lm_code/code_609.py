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
    Detects if the synthesis preserves a stereocenter throughout the route.

    This function checks if at least one stereocenter is preserved from the final product
    back through the synthesis route. It ignores starting materials (in_stock=True).
    """
    # Track if we found a preserved stereocenter
    preserved_stereocenter = False

    # Get the final product (root node)
    if route["type"] != "mol":
        print("Root node is not a molecule")
        return False

    # Get stereocenters in the final product
    try:
        final_mol = Chem.MolFromSmiles(route["smiles"])
        if final_mol is None:
            print(f"Could not parse final product SMILES: {route['smiles']}")
            return False

        final_stereocenters = Chem.FindMolChiralCenters(final_mol, includeUnassigned=False)

        if not final_stereocenters:
            print(f"Final product has no stereocenters: {route['smiles']}")
            return False

        print(f"Final product stereocenters: {final_stereocenters}")

        # Track stereocenters through the synthesis
        preserved_stereocenter = track_stereocenters(route, final_stereocenters)

    except Exception as e:
        print(f"Error analyzing final product: {e}")
        return False

    if preserved_stereocenter:
        print("At least one stereocenter is preserved throughout the synthesis")
    else:
        print("No stereocenter is preserved throughout the synthesis")

    return preserved_stereocenter

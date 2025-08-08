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


def main(node, target_stereocenters):
    """
    Recursively track stereocenters through the synthesis route.
    Returns True if at least one stereocenter is preserved throughout.
    """
    # Base case: starting material
    if node.get("in_stock", False):
        return False  # We don't expect starting materials to preserve stereocenters

    # If this is a molecule node, check for stereocenters
    if node["type"] == "mol":
        try:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol is None:
                print(f"Could not parse molecule SMILES: {node['smiles']}")
                return False

            stereocenters = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
            print(f"Molecule {node['smiles']} stereocenters: {stereocenters}")

            # If no stereocenters, this branch doesn't preserve stereocenters
            if not stereocenters:
                return False

            # Check if any of the target stereocenters are present in this molecule
            # This is a simplified check - in a real implementation, you would need to
            # track atom mappings through reactions to follow specific stereocenters
            for center, chirality in target_stereocenters:
                for mol_center, mol_chirality in stereocenters:
                    if mol_chirality == chirality:  # Simple check for same chirality
                        # For a more accurate check, we would need to track atom mappings
                        return True

            # If we have children, check if any branch preserves stereocenters
            for child in node.get("children", []):
                if track_stereocenters(child, stereocenters):
                    return True

            return False

        except Exception as e:
            print(f"Error analyzing molecule {node.get('smiles', 'unknown')}: {e}")
            return False

    # If this is a reaction node, check if any reactant preserves stereocenters
    elif node["type"] == "reaction":
        # For reaction nodes, check if any child (reactant) preserves stereocenters
        for child in node.get("children", []):
            if track_stereocenters(child, target_stereocenters):
                return True

        return False

    return False

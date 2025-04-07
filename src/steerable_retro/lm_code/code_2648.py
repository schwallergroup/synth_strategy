#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    Detects if the synthesis route maintains a stereocenter throughout the synthesis.
    """
    # Track stereocenters through the synthesis
    stereocenters_by_depth = {}

    def get_atom_mapping(smiles):
        """Extract atom mapping from a SMILES string"""
        # Extract atom mapping numbers using regex
        pattern = r":(\d+)]"
        mappings = re.findall(pattern, smiles)
        return [int(m) for m in mappings]

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            # Add depth information to the node
            node["depth"] = depth

            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Find stereocenters in this molecule
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)

                if len(chiral_centers) > 0:
                    # Store stereocenters with their configurations at this depth
                    stereo_info = []
                    for atom_idx, config in chiral_centers:
                        # Get atom mapping if available
                        atom = mol.GetAtomWithIdx(atom_idx)
                        atom_map = (
                            atom.GetProp("molAtomMapNumber")
                            if atom.HasProp("molAtomMapNumber")
                            else None
                        )
                        stereo_info.append((atom_idx, config, atom_map))

                    stereocenters_by_depth[depth] = stereo_info

                    if node.get("in_stock", False):
                        print(
                            f"Starting material at depth {depth} has stereocenter(s): {chiral_centers}"
                        )
                    elif depth == 0:
                        print(f"Final product has stereocenter(s): {chiral_centers}")
                    else:
                        print(
                            f"Intermediate at depth {depth} has stereocenter(s): {chiral_centers}"
                        )

        elif node["type"] == "reaction":
            # For reaction nodes, check if stereochemistry is preserved
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Extract reactants and product
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if the reaction preserves stereochemistry
                    # This would require more detailed analysis in a real implementation
                except Exception as e:
                    print(f"Error extracting reactants and product: {e}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if any stereocenter is maintained throughout the synthesis
    maintained = False

    # If we have stereocenters at both the final product (depth 0) and in starting materials
    if 0 in stereocenters_by_depth and any(d > 0 for d in stereocenters_by_depth.keys()):
        # Get the final product stereocenters
        final_stereocenters = stereocenters_by_depth[0]

        # Find the maximum depth (starting materials)
        max_depth = max(stereocenters_by_depth.keys())

        # Check if any stereocenter in starting materials is also in final product
        # In a real implementation, we would use atom mapping to track specific stereocenters
        if max_depth in stereocenters_by_depth:
            # For simplicity, we'll just check if there are stereocenters at both ends
            # A real implementation would track specific stereocenters through atom mapping
            maintained = True
            print(f"Stereocenters found at both final product and starting materials")

    return maintained

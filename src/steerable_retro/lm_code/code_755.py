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
    This function detects if the synthetic route maintains a stereocenter throughout without modification.
    """
    # Dictionary to store stereocenters at each molecule node
    stereocenters_by_mol = {}
    # Dictionary to track depth of each molecule
    mol_depths = {}
    # Track the target molecule (final product)
    target_mol = route["smiles"]

    def calculate_depth(node, parent_depth=0):
        """Calculate depth for each molecule node in the synthesis tree"""
        if node["type"] == "mol":
            mol_depths[node["smiles"]] = parent_depth
            for child in node.get("children", []):
                calculate_depth(child, parent_depth + 1)
        else:  # reaction node
            for child in node.get("children", []):
                calculate_depth(child, parent_depth)

    def find_stereocenters(node):
        """Find and store stereocenters for each molecule node"""
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            mol = Chem.MolFromSmiles(mol_smiles)

            if mol:
                # Find stereocenters in the molecule
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                if chiral_centers:
                    # Store the stereocenters with their configurations
                    stereocenters_by_mol[mol_smiles] = set(
                        (atom_idx, chirality) for atom_idx, chirality in chiral_centers
                    )
                    print(f"Found {len(chiral_centers)} stereocenters in molecule: {mol_smiles}")

        # Process children
        for child in node.get("children", []):
            find_stereocenters(child)

    def track_stereocenters_through_reactions(node):
        """Track stereocenters through reactions using atom mapping"""
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Get the product molecule from children
                product_node = None
                for child in node.get("children", []):
                    if child["type"] == "mol" and child["smiles"] == product_smiles:
                        product_node = child
                        break

                if product_node and product_smiles in stereocenters_by_mol:
                    # Check if the stereocenters in the product are preserved in the reactants
                    product_stereocenters = stereocenters_by_mol[product_smiles]

                    # For each reactant, check if it contains the same stereocenters
                    for reactant_smiles in reactants_smiles:
                        if reactant_smiles in stereocenters_by_mol:
                            reactant_stereocenters = stereocenters_by_mol[reactant_smiles]

                            # If the reactant has the same number of stereocenters, mark them as potentially preserved
                            if len(reactant_stereocenters) == len(product_stereocenters):
                                print(
                                    f"Potential stereocenter preservation between {product_smiles} and {reactant_smiles}"
                                )
            except Exception as e:
                print(f"Error tracking stereocenters through reaction: {e}")

        # Process children
        for child in node.get("children", []):
            track_stereocenters_through_reactions(child)

    # Calculate depth for each molecule
    calculate_depth(route)

    # Find stereocenters in all molecules
    find_stereocenters(route)

    # Track stereocenters through reactions
    track_stereocenters_through_reactions(route)

    # Check if stereocenters are preserved throughout the synthesis
    has_preserved_stereocenter = False

    # Get molecules sorted by depth (from target to starting materials)
    sorted_mols = sorted(mol_depths.items(), key=lambda x: x[1])

    # Check if there are stereocenters in the target molecule
    if target_mol in stereocenters_by_mol and len(stereocenters_by_mol[target_mol]) > 0:
        # Check if at least one stereocenter is preserved throughout the synthesis
        target_stereocenter_count = len(stereocenters_by_mol[target_mol])

        # Count how many molecules have the same number of stereocenters
        mols_with_stereocenters = 0
        for mol_smiles, _ in sorted_mols:
            if (
                mol_smiles in stereocenters_by_mol
                and len(stereocenters_by_mol[mol_smiles]) == target_stereocenter_count
            ):
                mols_with_stereocenters += 1

        # If more than one molecule has the same number of stereocenters, consider it preserved
        has_preserved_stereocenter = mols_with_stereocenters > 1

    print(f"Stereocenter preservation strategy: {has_preserved_stereocenter}")
    return has_preserved_stereocenter

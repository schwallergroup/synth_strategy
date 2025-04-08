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
    This function detects a synthesis strategy involving intramolecular ring formation
    creating a bicyclic system from an acyl chloride intermediate.
    """
    ring_formation_detected = False

    def dfs_traverse(node):
        nonlocal ring_formation_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Try to convert to RDKit molecules
                try:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if all(reactant_mols) and product_mol:
                        # Count rings in reactants and product
                        reactant_ring_count = sum(
                            len(mol.GetRingInfo().AtomRings()) for mol in reactant_mols
                        )
                        product_ring_count = len(product_mol.GetRingInfo().AtomRings())

                        # Check if product has more rings than reactants
                        if product_ring_count > reactant_ring_count:
                            # Check if acyl chloride pattern exists in any reactant
                            acyl_chloride_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#17]")
                            has_acyl_chloride = any(
                                mol.HasSubstructMatch(acyl_chloride_pattern)
                                for mol in reactant_mols
                            )

                            # Check if product has a bicyclic system
                            ring_systems = product_mol.GetRingInfo().AtomRings()
                            atoms_in_rings = set()
                            for ring in ring_systems:
                                atoms_in_rings.update(ring)

                            # Count atoms that appear in multiple rings
                            shared_atoms = 0
                            for atom_idx in atoms_in_rings:
                                rings_containing_atom = 0
                                for ring in ring_systems:
                                    if atom_idx in ring:
                                        rings_containing_atom += 1
                                if rings_containing_atom > 1:
                                    shared_atoms += 1

                            # If we have shared atoms between rings and had an acyl chloride, it's likely
                            # an intramolecular ring formation creating a bicyclic system
                            if shared_atoms > 0 and has_acyl_chloride:
                                print(
                                    f"Detected intramolecular ring formation creating a bicyclic system"
                                )
                                ring_formation_detected = True
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return ring_formation_detected

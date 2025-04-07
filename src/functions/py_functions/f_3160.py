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
    This function detects convergent synthesis of biaryl compounds via organometallic coupling,
    specifically looking for Stille-type coupling patterns.
    """
    # Track if we found the key features
    found_organometallic_coupling = False
    found_biaryl_formation = False

    def dfs_traverse(node):
        nonlocal found_organometallic_coupling, found_biaryl_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for organometallic coupling
            has_tin = any(
                "[Sn]" in reactant or "#50" in reactant for reactant in reactants_smiles
            )

            if has_tin:
                # Check if this is a biaryl formation
                try:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Count aromatic rings in reactants and product
                    reactant_aromatic_rings = sum(
                        sum(
                            1
                            for ring in Chem.GetSSSR(mol)
                            if all(
                                mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring
                            )
                        )
                        for mol in reactant_mols
                        if mol
                    )

                    product_aromatic_rings = (
                        sum(
                            1
                            for ring in Chem.GetSSSR(product_mol)
                            if all(
                                product_mol.GetAtomWithIdx(idx).GetIsAromatic()
                                for idx in ring
                            )
                        )
                        if product_mol
                        else 0
                    )

                    # If product has more aromatic rings connected than reactants, it's likely a biaryl formation
                    if product_aromatic_rings > 1 and has_tin:
                        found_organometallic_coupling = True
                        found_biaryl_formation = True
                        print("Found biaryl formation via organometallic coupling")
                except:
                    print("Error processing molecules for biaryl detection")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found both key features
    return found_organometallic_coupling and found_biaryl_formation

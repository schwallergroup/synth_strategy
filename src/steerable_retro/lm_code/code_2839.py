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
    This function detects the Gabriel synthesis pathway for amine formation:
    1. Activation of alcohol as mesylate
    2. Displacement with phthalimide
    3. Phthalimide deprotection to primary amine
    """
    # Initialize tracking variables
    has_alcohol_to_mesylate = False
    has_mesylate_to_phthalimide = False
    has_phthalimide_deprotection = False

    def dfs_traverse(node):
        nonlocal has_alcohol_to_mesylate, has_mesylate_to_phthalimide, has_phthalimide_deprotection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol to mesylate conversion
            alcohol_pattern = Chem.MolFromSmarts("[#8H1]")
            mesylate_pattern = Chem.MolFromSmarts("[#8][S](=[#8])(=[#8])[#6]")

            # Check for mesylate to phthalimide conversion
            phthalimide_pattern = Chem.MolFromSmarts(
                "[#8]=[#6]1[#6]2[#6][#6][#6][#6][#6]2[#6](=[#8])[#7]1"
            )

            # Check for phthalimide deprotection to primary amine
            primary_amine_pattern = Chem.MolFromSmarts("[#7H2]")

            # Convert SMILES to RDKit molecules
            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

                # Check for alcohol to mesylate
                if (
                    any(m and m.HasSubstructMatch(alcohol_pattern) for m in reactant_mols)
                    and product_mol
                    and product_mol.HasSubstructMatch(mesylate_pattern)
                ):
                    print("Detected alcohol to mesylate conversion")
                    has_alcohol_to_mesylate = True

                # Check for mesylate to phthalimide
                if (
                    any(m and m.HasSubstructMatch(mesylate_pattern) for m in reactant_mols)
                    and product_mol
                    and product_mol.HasSubstructMatch(phthalimide_pattern)
                ):
                    print("Detected mesylate to phthalimide conversion")
                    has_mesylate_to_phthalimide = True

                # Check for phthalimide deprotection
                if (
                    any(m and m.HasSubstructMatch(phthalimide_pattern) for m in reactant_mols)
                    and product_mol
                    and product_mol.HasSubstructMatch(primary_amine_pattern)
                ):
                    print("Detected phthalimide deprotection to primary amine")
                    has_phthalimide_deprotection = True

            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if all three steps of Gabriel synthesis are detected
    gabriel_synthesis_detected = (
        has_alcohol_to_mesylate and has_mesylate_to_phthalimide and has_phthalimide_deprotection
    )
    print(f"Gabriel synthesis pathway detected: {gabriel_synthesis_detected}")
    return gabriel_synthesis_detected

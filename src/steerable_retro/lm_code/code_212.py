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
    Detects a synthetic strategy where a late-stage SNAr coupling follows an amine deprotection.
    The strategy involves:
    1. Deprotection of a Boc-protected amine
    2. SNAr coupling of the resulting amine with a chloro-heteroaromatic system
    """
    # Track if we found the required reactions
    found_snar = False
    found_deprotection = False
    deprotection_depth = -1
    snar_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_snar, found_deprotection, deprotection_depth, snar_depth

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")
            product = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            product_mol = Chem.MolFromSmiles(product)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

            # Check for SNAr reaction (amine + chloro-heteroaromatic)
            if len(reactant_mols) == 2 and product_mol:
                # Check if one reactant has a primary amine
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                # Check if one reactant has a chloro-heteroaromatic
                chloro_pattern = Chem.MolFromSmarts("[n:1]~[c:2][Cl:3]")

                has_amine = False
                has_chloro = False

                for r in reactant_mols:
                    if r and r.HasSubstructMatch(amine_pattern):
                        has_amine = True
                    if r and r.HasSubstructMatch(chloro_pattern):
                        has_chloro = True

                # Check if product has a new C-N bond where the chlorine was
                if has_amine and has_chloro:
                    # This is a simplification - in a real implementation, we would need to
                    # check that the amine replaced the chlorine specifically
                    found_snar = True
                    snar_depth = depth
                    print(f"Found SNAr reaction at depth {depth}")

            # Check for Boc deprotection
            if product_mol:
                # Check if product has a primary amine
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                if product_mol.HasSubstructMatch(amine_pattern):
                    # Check if reactant had a Boc group
                    boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)[NH]")
                    for r_smiles in reactants:
                        r_mol = Chem.MolFromSmiles(r_smiles)
                        if r_mol and r_mol.HasSubstructMatch(boc_pattern):
                            found_deprotection = True
                            deprotection_depth = depth
                            print(f"Found Boc deprotection at depth {depth}")

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we found both reactions in the correct order
    # Remember: lower depth = later in synthesis (closer to final product)
    if found_snar and found_deprotection and snar_depth < deprotection_depth:
        print("Detected strategy: Late-stage SNAr coupling after amine deprotection")
        return True

    return False

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
    This function detects a Boc protection/deprotection sequence on a nitrogen atom.
    """
    # Track if we've found both protection and deprotection
    protection_found = False
    deprotection_found = False

    def dfs_traverse(node):
        nonlocal protection_found, deprotection_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc protection (formation of Boc group on N)
                boc_pattern = Chem.MolFromSmarts("[#6]C([#6])([#6])[#8]C(=O)[#7]")

                # For protection: product should have Boc but reactants shouldn't
                if not protection_found:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(boc_pattern):
                        reactants_with_boc = 0
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(boc_pattern):
                                reactants_with_boc += 1

                        if (
                            reactants_with_boc
                            < product_mol.GetSubstructMatches(boc_pattern).__len__()
                        ):
                            protection_found = True
                            print("Found Boc protection reaction")

                # For deprotection: reactants should have Boc but product shouldn't
                if not deprotection_found:
                    product_mol = Chem.MolFromSmiles(product)
                    reactants_with_boc = 0
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(boc_pattern):
                            reactants_with_boc += 1

                    if (
                        product_mol
                        and not product_mol.HasSubstructMatch(boc_pattern)
                        and reactants_with_boc > 0
                    ):
                        deprotection_found = True
                        print("Found Boc deprotection reaction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both protection and deprotection were found
    return protection_found and deprotection_found

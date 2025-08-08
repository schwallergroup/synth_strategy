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
    Detects nucleophilic aromatic substitution (SNAr) to form aryl ether linkage
    from aryl chloride and phenol.
    """
    snar_found = False

    def dfs_traverse(node):
        nonlocal snar_found

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            if len(reactants) >= 2:  # SNAr requires at least two reactants
                # Check for aryl chloride and phenol reactants
                aryl_chloride_pattern = Chem.MolFromSmarts("c-Cl")
                phenol_pattern = Chem.MolFromSmarts("c-[OH]")

                aryl_chloride_found = False
                phenol_found = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(aryl_chloride_pattern):
                            aryl_chloride_found = True
                        if reactant_mol.HasSubstructMatch(phenol_pattern):
                            phenol_found = True

                if aryl_chloride_found and phenol_found:
                    # Check if product has aryl ether linkage
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        aryl_ether_pattern = Chem.MolFromSmarts("c-O-c")
                        if product_mol.HasSubstructMatch(aryl_ether_pattern):
                            snar_found = True
                            print("Found SNAr reaction forming aryl ether linkage")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return snar_found

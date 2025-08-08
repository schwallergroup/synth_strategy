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
    This function detects a synthetic strategy involving nucleophilic aromatic substitution
    of a chlorine with an amine.
    """
    snAr_detected = False

    def dfs_traverse(node):
        nonlocal snAr_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for chloro-aromatic in reactants
                chloro_aromatic_pattern = Chem.MolFromSmarts("c[Cl]")
                chloro_present = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(chloro_aromatic_pattern):
                            chloro_present = True
                            print("Chloro-aromatic detected in reactants")
                    except:
                        continue

                # Check for amine reactant
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                amine_present = False
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(amine_pattern):
                            amine_present = True
                            print("Amine detected in reactants")
                    except:
                        continue

                # Check for C-N bond in product
                cn_bond_pattern = Chem.MolFromSmarts("c[NH]")
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(cn_bond_pattern):
                        if chloro_present and amine_present:
                            snAr_detected = True
                            print("Nucleophilic aromatic substitution detected")
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return snAr_detected

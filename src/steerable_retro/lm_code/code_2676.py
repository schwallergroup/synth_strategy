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
    Detects if the synthesis includes an SNAr reaction with an amine nucleophile.
    """
    has_snar_with_amine = False

    def dfs_traverse(node):
        nonlocal has_snar_with_amine

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Look for patterns indicating SNAr with amine
            # Check if one reactant has Cl attached to aromatic and another has NH2 or NH group
            has_aryl_chloride = False
            has_amine = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    aryl_cl_pattern = Chem.MolFromSmarts("c[Cl]")
                    amine_pattern = Chem.MolFromSmarts("[NH2,NH]")

                    if reactant_mol.HasSubstructMatch(aryl_cl_pattern):
                        has_aryl_chloride = True
                    if reactant_mol.HasSubstructMatch(amine_pattern):
                        has_amine = True

            # Check if product has C-N bond where Cl was
            if has_aryl_chloride and has_amine:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    c_n_pattern = Chem.MolFromSmarts("c[N]")
                    if product_mol.HasSubstructMatch(c_n_pattern):
                        has_snar_with_amine = True
                        print("Detected SNAr with amine")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_snar_with_amine

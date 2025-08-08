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
    This function detects if the synthetic route involves reductive amination to connect fragments.
    """
    reductive_amination_detected = False

    def dfs_traverse(node):
        nonlocal reductive_amination_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for patterns consistent with reductive amination
                # 1. One reactant has C=O (ketone/aldehyde)
                # 2. Another reactant has NH or NH2
                # 3. Product has new C-N bond where C=O was

                ketone_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#6]")
                amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")

                has_ketone = False
                has_amine = False

                for reactant_smiles in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(ketone_pattern):
                            has_ketone = True
                        if reactant_mol.HasSubstructMatch(amine_pattern):
                            has_amine = True

                if has_ketone and has_amine:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol:
                        # Check if product has a new C-N bond
                        cn_bond_pattern = Chem.MolFromSmarts("[#6]-[#7]")
                        if product_mol.HasSubstructMatch(cn_bond_pattern):
                            reductive_amination_detected = True
                            print("Reductive amination detected in reaction:", rsmi)

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return reductive_amination_detected

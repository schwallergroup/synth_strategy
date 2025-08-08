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
    Detects if the synthesis route involves a late-stage fragment coupling via S-alkylation,
    where a thiol (SH) reacts with a chloromethyl group to form a thioether linkage.
    """
    s_alkylation_found = False
    final_step = True

    def dfs_traverse(node):
        nonlocal s_alkylation_found, final_step

        if node["type"] == "reaction" and final_step:
            # This is the first reaction we encounter (latest stage)
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if one reactant has a thiol group
            thiol_pattern = Chem.MolFromSmarts("[#16H1]")
            chloromethyl_pattern = Chem.MolFromSmarts("[Cl][CH2][c]")
            thioether_pattern = Chem.MolFromSmarts("[#16][CH2][c]")

            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol and product_mol.HasSubstructMatch(thioether_pattern):
                thiol_found = False
                chloromethyl_found = False

                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(thiol_pattern):
                            thiol_found = True
                        if reactant_mol.HasSubstructMatch(chloromethyl_pattern):
                            chloromethyl_found = True

                if thiol_found and chloromethyl_found:
                    print("Found late-stage S-alkylation coupling")
                    s_alkylation_found = True

            final_step = False

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return s_alkylation_found

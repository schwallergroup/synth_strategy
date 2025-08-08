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
    This function detects if the synthetic route involves a late-stage mesylation
    (conversion of alcohol to methanesulfonate) as the final step.
    """
    mesylation_detected = False
    final_step = True

    def dfs_traverse(node):
        nonlocal mesylation_detected, final_step

        if node["type"] == "reaction" and final_step:
            # Get reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Create RDKit mol objects
            product_mol = Chem.MolFromSmiles(product)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

            # Check if any reactant has an OH group
            alcohol_reactant = False
            for r_mol in reactant_mols:
                if r_mol and r_mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]")):
                    alcohol_reactant = True
                    break

            # Check if product has a mesylate group
            mesylate_pattern = Chem.MolFromSmarts("[O][S](=[O])(=[O])[C]")
            if product_mol and alcohol_reactant and product_mol.HasSubstructMatch(mesylate_pattern):
                mesylation_detected = True

            final_step = False

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    print(f"Late-stage mesylation strategy detected: {mesylation_detected}")
    return mesylation_detected

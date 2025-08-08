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
    This function detects a synthetic strategy involving reductive amination.
    """
    reductive_amination_detected = False

    def dfs_traverse(node):
        nonlocal reductive_amination_detected

        if node["type"] == "reaction" and not reductive_amination_detected:
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check for reductive amination pattern:
                    # 1. One reactant has C=O (ketone/aldehyde)
                    # 2. One reactant has NH or NH2
                    # 3. Product has C-N where C was part of C=O

                    carbonyl_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#6]")  # Ketone
                    amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]")  # Primary or secondary amine

                    product_mol = Chem.MolFromSmiles(product)

                    has_carbonyl = False
                    has_amine = False

                    for r in reactants:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol:
                            if r_mol.HasSubstructMatch(carbonyl_pattern):
                                has_carbonyl = True
                            if r_mol.HasSubstructMatch(amine_pattern):
                                has_amine = True

                    # If product has tertiary amine and reactants had carbonyl and amine
                    if has_carbonyl and has_amine:
                        tertiary_amine_pattern = Chem.MolFromSmarts("[#6]-[#7](-[#6])-[#6]")
                        if product_mol and product_mol.HasSubstructMatch(tertiary_amine_pattern):
                            reductive_amination_detected = True
                            print(f"Reductive amination detected in reaction: {rsmi}")
                except:
                    pass

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return reductive_amination_detected

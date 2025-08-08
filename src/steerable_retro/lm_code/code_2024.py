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
    This function detects if the synthesis includes multiple amide bond formation steps.
    """
    amide_formation_count = 0

    def dfs_traverse(node):
        nonlocal amide_formation_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carboxylic acid in reactants
                acid_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#8;H1,-]")
                amine_pattern = Chem.MolFromSmarts("[#7;!$(NC=O);!$(N=*);!$([#7]-[#7])]")
                amide_pattern = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#7]")

                # Check if reactants contain acid and amine, and product contains amide
                reactants_have_acid = any(
                    Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).HasSubstructMatch(acid_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                )
                reactants_have_amine = any(
                    Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).HasSubstructMatch(amine_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                )
                product_has_amide = Chem.MolFromSmiles(product) and Chem.MolFromSmiles(
                    product
                ).HasSubstructMatch(amide_pattern)

                if reactants_have_acid and reactants_have_amine and product_has_amide:
                    amide_formation_count += 1
                    print(f"Detected amide bond formation. Total count: {amide_formation_count}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return amide_formation_count >= 2  # Return True if at least 2 amide formations are detected

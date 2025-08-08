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
    Detects if the synthesis route involves a ketone reduction that creates a stereocenter.
    """
    has_ketone_reduction_with_stereocenter = False

    def dfs_traverse(node):
        nonlocal has_ketone_reduction_with_stereocenter

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check for ketone pattern in reactants
            ketone_pattern = Chem.MolFromSmarts("[C](=[O])[C]")

            # Check for secondary alcohol pattern in product
            alcohol_pattern = Chem.MolFromSmarts("[C]([OH])[C]")

            reactant = Chem.MolFromSmiles(reactants_smiles)
            product = Chem.MolFromSmiles(product_smiles)

            if reactant is not None and product is not None:
                ketone_present = reactant.HasSubstructMatch(ketone_pattern)
                alcohol_present = product.HasSubstructMatch(alcohol_pattern)

                # Check if product has a stereocenter
                has_stereocenter = False
                if alcohol_present:
                    for atom in product.GetAtoms():
                        if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                            has_stereocenter = True
                            break

                if ketone_present and alcohol_present and has_stereocenter:
                    has_ketone_reduction_with_stereocenter = True
                    print(
                        f"Detected ketone reduction with stereocenter formation in reaction: {node.get('metadata', {}).get('ID', '')}"
                    )

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_ketone_reduction_with_stereocenter

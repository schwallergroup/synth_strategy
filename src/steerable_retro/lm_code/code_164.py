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
    This function detects if the synthesis involves a trifluoromethoxy group in the final product.
    """
    trifluoromethoxy_pattern = Chem.MolFromSmarts("[#8][#6]([#9])([#9])[#9]")
    final_product_has_trifluoromethoxy = False

    if route["type"] == "mol":
        mol = Chem.MolFromSmiles(route["smiles"])
        if mol and mol.HasSubstructMatch(trifluoromethoxy_pattern):
            final_product_has_trifluoromethoxy = True

    print(f"Trifluoromethoxy in final product: {final_product_has_trifluoromethoxy}")
    return final_product_has_trifluoromethoxy

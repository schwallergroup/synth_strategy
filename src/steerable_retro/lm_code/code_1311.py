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
    Detects if the synthetic route involves an SNAr transformation (fluorine displacement).
    """
    snar_found = False

    def dfs_traverse(node):
        nonlocal snar_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for fluorinated aromatic in reactants
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        fluoro_aromatic = Chem.MolFromSmarts("c[F]")
                        if reactant_mol.HasSubstructMatch(fluoro_aromatic):
                            # Check for O-C bond formation in product
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol:
                                ether_pattern = Chem.MolFromSmarts("cOC")
                                if product_mol.HasSubstructMatch(ether_pattern):
                                    snar_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"SNAr transformation detected: {snar_found}")
    return snar_found

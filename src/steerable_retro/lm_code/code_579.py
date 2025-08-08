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
    This function detects a synthetic strategy involving C-C bond disconnection between
    an aromatic ring and a spiro carbon center.
    """
    has_cc_disconnection = False

    def dfs_traverse(node):
        nonlocal has_cc_disconnection

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check for C-C disconnection pattern
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                product = Chem.MolFromSmiles(product_smiles)

                if product and all(r for r in reactants):
                    # Pattern for C-C bond between aromatic and spiro carbon
                    cc_bond_pattern = Chem.MolFromSmarts("[c]-[C;R]")

                    # Check if product has the pattern but reactants don't have it connected
                    if product.HasSubstructMatch(cc_bond_pattern):
                        # Check if we have separate aromatic and aliphatic fragments in reactants
                        has_aromatic = False
                        has_aliphatic = False

                        for r in reactants:
                            if r:
                                aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
                                aliphatic_pattern = Chem.MolFromSmarts("[C;R][C;R][C;R][C;R]")

                                if r.HasSubstructMatch(aromatic_pattern):
                                    has_aromatic = True
                                if r.HasSubstructMatch(aliphatic_pattern):
                                    has_aliphatic = True

                        if has_aromatic and has_aliphatic:
                            has_cc_disconnection = True
                            print(
                                "Detected C-C bond disconnection between aromatic ring and spiro carbon"
                            )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_cc_disconnection

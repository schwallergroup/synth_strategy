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
    This function detects a synthetic strategy involving reduction of an aromatic
    N-heterocycle (like pyridine) to a saturated heterocycle (like piperidine).
    """
    has_heterocycle_reduction = False

    def dfs_traverse(node):
        nonlocal has_heterocycle_reduction

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for pyridine pattern in reactants
                pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")
                # Check for piperidine pattern in product
                piperidine_pattern = Chem.MolFromSmarts("C1CCNCC1")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if (
                    product_mol
                    and any(r and r.HasSubstructMatch(pyridine_pattern) for r in reactant_mols)
                    and product_mol.HasSubstructMatch(piperidine_pattern)
                ):
                    print("Detected heterocycle reduction: pyridine to piperidine")
                    has_heterocycle_reduction = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return has_heterocycle_reduction

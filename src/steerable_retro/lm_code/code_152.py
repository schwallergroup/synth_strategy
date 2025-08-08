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
    Detects the transformation of an amine to a nitrile group on an aryl ring.
    """
    # Track if we found the key patterns
    amine_to_nitrile_found = False

    # SMARTS patterns
    aryl_amine_pattern = Chem.MolFromSmarts("c-[NH2]")
    aryl_nitrile_pattern = Chem.MolFromSmarts("c-[C]#[N]")

    def dfs_traverse(node):
        nonlocal amine_to_nitrile_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactants = reactants_part.split(".")
                product = Chem.MolFromSmiles(product_part)

                # Check if any reactant has aryl-amine and product has aryl-nitrile
                if product and product.HasSubstructMatch(aryl_nitrile_pattern):
                    for r in reactants:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol and r_mol.HasSubstructMatch(aryl_amine_pattern):
                            amine_to_nitrile_found = True
                            print(f"Found amine to nitrile transformation: {r} -> {product_part}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    if amine_to_nitrile_found:
        print("Detected amine to nitrile transformation strategy")

    return amine_to_nitrile_found

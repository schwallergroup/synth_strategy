#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects if the synthetic route uses a protected amine/piperazine
    in an SNAr reaction.
    """
    found_protected_amine_snar = False

    def dfs_traverse(node, depth=0):
        nonlocal found_protected_amine_snar

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for SNAr pattern with protected amine
                fluoro_aromatic = Chem.MolFromSmarts("c([F])")
                protected_amine = Chem.MolFromSmarts("[#7]C(=[O])[O]C")  # Boc or similar protection

                # Check for C-N bond in product
                cn_bond_product = Chem.MolFromSmarts("c([#7])")

                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(cn_bond_product):
                    fluoro_found = False
                    protected_amine_found = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if not reactant_mol:
                            continue

                        if reactant_mol.HasSubstructMatch(fluoro_aromatic):
                            fluoro_found = True

                        if reactant_mol.HasSubstructMatch(protected_amine):
                            protected_amine_found = True

                    if fluoro_found and protected_amine_found:
                        print(f"Found protected amine in SNAr at depth {depth}")
                        found_protected_amine_snar = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return found_protected_amine_snar

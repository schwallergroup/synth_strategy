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
    This function detects if the synthetic route contains an SNAr reaction
    with a fluorinated aromatic compound.
    """
    found_snar = False

    def dfs_traverse(node, depth=0):
        nonlocal found_snar

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for SNAr pattern: fluorinated aromatic + nucleophile → substituted product
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    product_mol = Chem.MolFromSmiles(product)

                    if reactant_mol and product_mol:
                        # Pattern for fluorinated aromatic with electron-withdrawing group
                        fluoro_aromatic = Chem.MolFromSmarts("c([F])")
                        nitro_aromatic = Chem.MolFromSmarts("c([N+](=[O])[O-])")

                        # Pattern for nitrogen nucleophile in other reactant
                        nitrogen_nucleophile = Chem.MolFromSmarts("[#7;H,H2]")

                        # Pattern for C-N bond in product where F was
                        cn_bond_product = Chem.MolFromSmarts("c([#7])")

                        if (
                            reactant_mol.HasSubstructMatch(fluoro_aromatic)
                            and reactant_mol.HasSubstructMatch(nitro_aromatic)
                            and product_mol.HasSubstructMatch(cn_bond_product)
                        ):

                            # Check if any other reactant has nitrogen nucleophile
                            for other_reactant in reactants:
                                if other_reactant != reactant:
                                    other_mol = Chem.MolFromSmiles(other_reactant)
                                    if other_mol and other_mol.HasSubstructMatch(
                                        nitrogen_nucleophile
                                    ):
                                        print(
                                            f"Found SNAr with fluorinated aromatic at depth {depth}"
                                        )
                                        found_snar = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return found_snar

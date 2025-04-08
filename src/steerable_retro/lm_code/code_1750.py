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
    Detects if the synthetic route involves coupling of two halogenated fragments
    to form a biaryl amine.
    """

    def dfs_traverse(node):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for coupling of halogenated fragments
                fluoro_aromatic_pattern = Chem.MolFromSmarts("[c]-[F]")
                bromo_aromatic_pattern = Chem.MolFromSmarts("[c]-[Br]")
                chloro_aromatic_pattern = Chem.MolFromSmarts("[c]-[Cl]")
                biaryl_amine_pattern = Chem.MolFromSmarts("[c]-[NH]-[c]")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                # Count halogenated fragments in reactants
                halogenated_fragments = 0
                for r in reactant_mols:
                    if r:
                        if (
                            r.HasSubstructMatch(fluoro_aromatic_pattern)
                            or r.HasSubstructMatch(bromo_aromatic_pattern)
                            or r.HasSubstructMatch(chloro_aromatic_pattern)
                        ):
                            halogenated_fragments += 1

                # Check if product forms biaryl amine
                forms_biaryl_amine = product_mol and product_mol.HasSubstructMatch(
                    biaryl_amine_pattern
                )

                if halogenated_fragments >= 2 and forms_biaryl_amine:
                    print(
                        f"Found coupling of {halogenated_fragments} halogenated fragments to form biaryl amine"
                    )
                    return True

        # Traverse children
        for child in node.get("children", []):
            if dfs_traverse(child):
                return True

        return False

    # Start traversal from the root
    return dfs_traverse(route)

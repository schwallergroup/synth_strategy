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
    Detects if the synthetic route involves a late-stage carbamoylation (depth 0-1).
    """
    late_carbamoylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_carbamoylation_found

        if node["type"] == "reaction" and depth <= 1:  # Only check at depths 0-1 (late stage)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for secondary amine in reactants
                secondary_amine_in_reactants = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        secondary_amine = Chem.MolFromSmarts("[#7;H1]")
                        if reactant_mol.HasSubstructMatch(secondary_amine):
                            secondary_amine_in_reactants = True

                # Check for carbamate in product
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and secondary_amine_in_reactants:
                    carbamate = Chem.MolFromSmarts("[#7]C(=[#8])[#8]C")
                    if product_mol.HasSubstructMatch(carbamate):
                        late_carbamoylation_found = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Late-stage carbamoylation detected: {late_carbamoylation_found}")
    return late_carbamoylation_found

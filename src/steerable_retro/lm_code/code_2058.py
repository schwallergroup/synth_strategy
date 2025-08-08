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
    Detects a nitro reduction to amine in the synthesis route.
    """
    has_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal has_nitro_reduction

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro group in reactants
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")

                    if (
                        reactant_mol
                        and nitro_pattern
                        and reactant_mol.HasSubstructMatch(nitro_pattern)
                    ):
                        # Check for amine in product
                        product_mol = Chem.MolFromSmiles(product)
                        amine_pattern = Chem.MolFromSmarts("[NH2]")

                        if (
                            product_mol
                            and amine_pattern
                            and product_mol.HasSubstructMatch(amine_pattern)
                        ):
                            has_nitro_reduction = True
                            print(
                                f"Detected nitro reduction at depth {node.get('depth', 'unknown')}"
                            )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_nitro_reduction

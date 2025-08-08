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
    This function detects a strategy involving sequential aromatic functionalization
    (bromination followed by coupling or nitro reduction).
    """
    aromatic_bromination = False
    nitro_reduction = False

    def dfs_traverse(node):
        nonlocal aromatic_bromination, nitro_reduction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            try:
                # Parse reactants and products
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                product = Chem.MolFromSmiles(product_smiles)

                if all(r is not None for r in reactants) and product is not None:
                    # Check for aromatic bromination
                    aryl_br_pattern = Chem.MolFromSmarts("[#6;a]-[Br]")

                    product_has_aryl_br = product.HasSubstructMatch(aryl_br_pattern)
                    reactants_have_aryl_br = any(
                        r.HasSubstructMatch(aryl_br_pattern) for r in reactants if r is not None
                    )

                    if product_has_aryl_br and not reactants_have_aryl_br:
                        aromatic_bromination = True
                        print(f"Aromatic bromination detected at reaction with RSMI: {rsmi}")

                    # Check for nitro reduction
                    nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])-[O-]")
                    amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")

                    reactants_have_nitro = any(
                        r.HasSubstructMatch(nitro_pattern) for r in reactants if r is not None
                    )
                    product_has_amine = product.HasSubstructMatch(amine_pattern)

                    if reactants_have_nitro and product_has_amine:
                        nitro_reduction = True
                        print(f"Nitro reduction detected at reaction with RSMI: {rsmi}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if the strategy is detected
    return aromatic_bromination and nitro_reduction

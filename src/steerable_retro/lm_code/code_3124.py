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
    This function detects protection-deprotection cycling strategy with nitrogen protecting groups.
    It looks for multiple protection and deprotection events in the synthesis route.
    """
    protection_count = 0
    deprotection_count = 0

    # SMARTS patterns for protecting groups
    boc_pattern = Chem.MolFromSmarts("[#6](C)(C)(C)OC(=O)[#7]")
    acetyl_pattern = Chem.MolFromSmarts("CC(=O)[#7]")

    def dfs_traverse(node):
        nonlocal protection_count, deprotection_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            products_smiles = rsmi.split(">")[-1]

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".") if r]
                product = Chem.MolFromSmiles(products_smiles)

                # Check for protection (appearance of protecting group)
                reactant_has_boc = any(
                    r is not None and r.HasSubstructMatch(boc_pattern) for r in reactants
                )
                reactant_has_acetyl = any(
                    r is not None and r.HasSubstructMatch(acetyl_pattern) for r in reactants
                )
                product_has_boc = product is not None and product.HasSubstructMatch(boc_pattern)
                product_has_acetyl = product is not None and product.HasSubstructMatch(
                    acetyl_pattern
                )

                if (not reactant_has_boc and product_has_boc) or (
                    not reactant_has_acetyl and product_has_acetyl
                ):
                    protection_count += 1
                    print(f"Protection detected: {rsmi}")

                # Check for deprotection (disappearance of protecting group)
                if (reactant_has_boc and not product_has_boc) or (
                    reactant_has_acetyl and not product_has_acetyl
                ):
                    deprotection_count += 1
                    print(f"Deprotection detected: {rsmi}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we have at least one protection and one deprotection
    has_cycles = protection_count >= 1 and deprotection_count >= 1
    print(f"Protection count: {protection_count}, Deprotection count: {deprotection_count}")
    return has_cycles

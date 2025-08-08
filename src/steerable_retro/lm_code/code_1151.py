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
    Detects if the synthetic route contains protection-deprotection sequences,
    particularly focusing on Boc protection/deprotection and nitro reduction.
    """
    has_boc_protection = False
    has_boc_deprotection = False
    has_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal has_boc_protection, has_boc_deprotection, has_nitro_reduction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc protection
            if any("C(=O)OC(C)(C)C" not in r for r in reactants) and "C(=O)OC(C)(C)C" in product:
                print("Found Boc protection")
                has_boc_protection = True

            # Check for Boc deprotection
            if "C(=O)OC(C)(C)C" in "".join(reactants) and "C(=O)OC(C)(C)C" not in product:
                print("Found Boc deprotection")
                has_boc_deprotection = True

            # Check for nitro reduction
            nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
            amine_pattern = Chem.MolFromSmarts("[NX3;H2]")

            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if (
                    any(mol and mol.HasSubstructMatch(nitro_pattern) for mol in reactant_mols)
                    and product_mol
                    and product_mol.HasSubstructMatch(amine_pattern)
                    and not product_mol.HasSubstructMatch(nitro_pattern)
                ):
                    print("Found nitro reduction")
                    has_nitro_reduction = True
            except:
                pass  # Handle parsing errors gracefully

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if at least two protection/deprotection operations are found
    return (
        (has_boc_protection and has_boc_deprotection)
        or has_nitro_reduction
        or (has_boc_protection and has_nitro_reduction)
        or (has_boc_deprotection and has_nitro_reduction)
    )

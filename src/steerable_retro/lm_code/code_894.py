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
    Detects a synthetic strategy involving protection/deprotection sequences.
    Focuses on N-acetylation and tetrahydropyran protecting groups.
    """
    acetyl_protection = False
    thp_protection = False

    def dfs_traverse(node, depth=0):
        nonlocal acetyl_protection, thp_protection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".") if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and all(r for r in reactants):
                # Check for N-acetylation
                acetyl_pattern = Chem.MolFromSmarts("[N]-[C](=[O])-[CH3]")

                product_has_acetyl = product.HasSubstructMatch(acetyl_pattern)
                reactants_have_acetyl = any(r.HasSubstructMatch(acetyl_pattern) for r in reactants)

                if (product_has_acetyl and not reactants_have_acetyl) or (
                    not product_has_acetyl and reactants_have_acetyl
                ):
                    acetyl_protection = True
                    print(f"N-acetyl protection/deprotection detected at depth {depth}")

                # Check for tetrahydropyran protection
                thp_pattern = Chem.MolFromSmarts("[CH2]1[CH2][CH2][CH2][CH2]O1")

                product_has_thp = product.HasSubstructMatch(thp_pattern)
                reactants_have_thp = any(r.HasSubstructMatch(thp_pattern) for r in reactants)

                if (product_has_thp and not reactants_have_thp) or (
                    not product_has_thp and reactants_have_thp
                ):
                    thp_protection = True
                    print(f"Tetrahydropyran protection/deprotection detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = acetyl_protection or thp_protection
    print(f"Protection/deprotection sequence: {result}")
    return result

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
    Detects if the synthetic route involves protection of carboxylic acid as
    methyl ester followed by deprotection later in the synthesis.
    """
    # Track if we've seen protection and deprotection
    protection_found = False
    deprotection_found = False

    def dfs_traverse(node):
        nonlocal protection_found, deprotection_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for esterification (protection)
                carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")
                methyl_ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product) if product else None

                if (
                    product_mol
                    and any(
                        r and r.HasSubstructMatch(carboxylic_acid_pattern) for r in reactant_mols
                    )
                    and product_mol.HasSubstructMatch(methyl_ester_pattern)
                ):
                    protection_found = True
                    print("Found carboxylic acid protection to methyl ester")

                # Check for hydrolysis (deprotection)
                if (
                    any(r and r.HasSubstructMatch(methyl_ester_pattern) for r in reactant_mols)
                    and product_mol
                    and product_mol.HasSubstructMatch(carboxylic_acid_pattern)
                ):
                    deprotection_found = True
                    print("Found methyl ester deprotection to carboxylic acid")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both protection and deprotection were found
    return protection_found and deprotection_found

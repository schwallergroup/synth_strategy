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
    This function detects if the synthetic route involves formation of an ether linkage.
    """
    ether_formation_detected = False

    def dfs_traverse(node):
        nonlocal ether_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for alcohol in reactants and ether in product
                alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")
                ether_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol and all(r for r in reactant_mols):
                    has_alcohol = any(
                        r.HasSubstructMatch(alcohol_pattern) for r in reactant_mols if r
                    )
                    product_has_ether = product_mol.HasSubstructMatch(ether_pattern)

                    if has_alcohol and product_has_ether:
                        # Count ethers in reactants vs product to confirm formation
                        reactant_ether_count = sum(
                            len(r.GetSubstructMatches(ether_pattern)) for r in reactant_mols if r
                        )
                        product_ether_count = len(product_mol.GetSubstructMatches(ether_pattern))

                        if product_ether_count > reactant_ether_count:
                            print(f"Ether linkage formation detected in reaction: {rsmi}")
                            ether_formation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return ether_formation_detected

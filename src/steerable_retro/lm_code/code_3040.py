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
    This function detects a strategy involving multiple benzylation reactions.
    It counts the number of benzylation reactions in the route.
    """
    benzylation_count = 0

    def dfs_traverse(node):
        nonlocal benzylation_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for benzylation: formation of benzyl ether
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol:
                # SMARTS for benzyl ether
                benzyl_ether_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]-c1ccccc1")

                if product_mol.HasSubstructMatch(benzyl_ether_pattern):
                    # Check if any reactant contains benzyl group attached to a leaving group
                    for r_smi in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r_smi)
                        if not r_mol:
                            continue

                        # Check for benzyl halide or similar
                        if r_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("c1ccccc1-[#6]-[F,Cl,Br,I,O]")
                        ):
                            print("Benzylation reaction detected:", rsmi)
                            benzylation_count += 1
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if there are at least 2 benzylation reactions
    return benzylation_count >= 2

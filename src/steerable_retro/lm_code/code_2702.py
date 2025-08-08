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
    Detects a synthetic strategy that uses Williamson ether synthesis for fragment coupling.
    """
    has_williamson_ether = False

    def dfs_traverse(node):
        nonlocal has_williamson_ether

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Williamson ether synthesis
                if len(reactants) >= 2:
                    benzyl_bromide_patt = Chem.MolFromSmarts("c1ccccc1-[#6]-[Br]")
                    alcohol_patt = Chem.MolFromSmarts("[#6]-[#8H]")
                    benzyl_ether_patt = Chem.MolFromSmarts("c1ccccc1-[#6]-[#8]-[#6]")

                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    product_mol = Chem.MolFromSmiles(product) if product else None

                    has_benzyl_bromide = any(
                        mol and mol.HasSubstructMatch(benzyl_bromide_patt) for mol in reactant_mols
                    )
                    has_alcohol = any(
                        mol and mol.HasSubstructMatch(alcohol_patt) for mol in reactant_mols
                    )

                    if (
                        has_benzyl_bromide
                        and has_alcohol
                        and product_mol
                        and product_mol.HasSubstructMatch(benzyl_ether_patt)
                    ):
                        print("Detected Williamson ether synthesis for fragment coupling")
                        has_williamson_ether = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_williamson_ether

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
    This function detects a strategy where a biaryl system is formed where one of the
    coupling partners contains a trifluoromethoxy group.
    """
    # Initialize tracking variables
    has_biaryl_coupling = False
    has_trifluoromethoxy_reactant = False

    def dfs_traverse(node, depth=0):
        nonlocal has_biaryl_coupling, has_trifluoromethoxy_reactant

        # Check for reactions
        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for biaryl coupling
            if any("B(O)" in r for r in reactants) and any("Br" in r for r in reactants):
                product_mol = Chem.MolFromSmiles(product)

                # Check for trifluoromethoxy in reactants
                for r in reactants:
                    reactant_mol = Chem.MolFromSmiles(r)
                    if reactant_mol:
                        trifluoromethoxy_pattern = Chem.MolFromSmarts("[O][C]([F])([F])[F]")
                        if reactant_mol.HasSubstructMatch(trifluoromethoxy_pattern):
                            has_trifluoromethoxy_reactant = True
                            print(f"Found trifluoromethoxy group in reactant: {r}")

                # Check if this is a biaryl formation
                if product_mol:
                    biaryl_pattern = Chem.MolFromSmarts(
                        "[c]!@[c]"
                    )  # Non-ring bond between two aromatic carbons
                    if product_mol.HasSubstructMatch(biaryl_pattern):
                        has_biaryl_coupling = True
                        print(f"Found biaryl coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = has_biaryl_coupling and has_trifluoromethoxy_reactant

    if strategy_present:
        print("Detected trifluoromethoxy-containing biaryl synthesis strategy")
    else:
        print("Strategy not detected")

    return strategy_present

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
    This function detects a convergent synthesis strategy with biaryl ether formation.
    It looks for a reaction where two complex fragments are joined via C-O-C bond formation.
    """
    has_biaryl_ether_formation = False

    def dfs_traverse(node):
        nonlocal has_biaryl_ether_formation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if we have multiple reactants (convergent)
                if len(reactants) >= 2:
                    # Check for biaryl ether formation
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # Look for biaryl ether pattern in product
                        biaryl_ether_pattern = Chem.MolFromSmarts("[c]-[O]-[c]")
                        if product_mol.HasSubstructMatch(biaryl_ether_pattern):
                            # Check if reactants have aromatic rings
                            aromatic_reactants = 0
                            for reactant in reactants:
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol and reactant_mol.HasSubstructMatch(
                                    Chem.MolFromSmarts("c")
                                ):
                                    aromatic_reactants += 1

                            # If at least 2 aromatic reactants, likely biaryl ether formation
                            if aromatic_reactants >= 2:
                                print("Found biaryl ether formation in convergent step")
                                has_biaryl_ether_formation = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_biaryl_ether_formation

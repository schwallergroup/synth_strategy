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
    Detects if the synthetic route involves construction of a heterocyclic scaffold,
    which is a common strategy in medicinal chemistry.
    """
    # Track ring formation reactions
    ring_formations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Count rings in reactants and product
                product_mol = Chem.MolFromSmiles(product)
                if not product_mol:
                    return

                product_rings = Chem.GetSSSR(product_mol)
                product_ring_count = len(product_rings)

                reactant_ring_count = 0
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        reactant_rings = Chem.GetSSSR(reactant_mol)
                        reactant_ring_count += len(reactant_rings)

                # If product has more rings than reactants combined, ring formation occurred
                if product_ring_count > reactant_ring_count:
                    # Check if the new ring contains heteroatoms
                    heteroatom_pattern = Chem.MolFromSmarts("[#7,#8,#16]")
                    if product_mol.HasSubstructMatch(heteroatom_pattern):
                        print(f"Found heterocycle formation at depth {depth}")
                        ring_formations.append(depth)
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # If we found at least one ring formation, the strategy is present
    return len(ring_formations) > 0

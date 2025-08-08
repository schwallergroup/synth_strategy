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
    This function detects early-stage formation of a nitrogen-containing heterocyclic ring.
    """
    has_ring_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_ring_formation

        if node["type"] == "reaction" and depth >= 3:  # Early stage (high depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for nitrogen-containing ring in product but not in reactants
                piperidine_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6][#7]1")

                has_ring_in_reactants = False
                for reactant in reactants_smiles:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(piperidine_pattern):
                            has_ring_in_reactants = True
                            break
                    except:
                        continue

                has_ring_in_product = False
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(piperidine_pattern):
                        has_ring_in_product = True
                except:
                    pass

                if not has_ring_in_reactants and has_ring_in_product:
                    print(f"Found heterocyclic ring formation at depth {depth}")
                    has_ring_formation = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_ring_formation

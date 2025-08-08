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
    This function detects the use of a formamidine reagent
    for heterocyclic ring formation, particularly in late-stage synthesis.
    """

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Check for formamidine-like structure (N=CH-NH2)
            # This is a simplification - in reality would need a more specific pattern
            formamidine_pattern = Chem.MolFromSmarts("[N]=[CH][NH2]")
            formamidine_in_reactants = any(
                r and len(r.GetSubstructMatches(formamidine_pattern)) > 0 for r in reactant_mols
            )

            if formamidine_in_reactants:
                # Check if product has a new ring
                reactants_ring_count = sum(
                    r.GetRingInfo().NumRings() if r else 0 for r in reactant_mols
                )
                product_ring_count = product_mol.GetRingInfo().NumRings() if product_mol else 0

                if product_ring_count > reactants_ring_count:
                    # Check if this is late-stage (depth <= 1)
                    if depth <= 1:
                        print(f"Late-stage formamidine-based cyclization detected at depth {depth}")
                        return True

        # Traverse children
        result = False
        for child in node.get("children", []):
            if dfs_traverse(child, depth + 1):
                result = True

        return result

    # Start traversal
    return dfs_traverse(route)

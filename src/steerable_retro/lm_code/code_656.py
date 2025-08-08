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
    This function detects if the synthesis involves a late-stage aromatic amine coupling.
    """
    late_stage_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_coupling

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if we have at least 2 reactants
            if len(reactants_smiles) >= 2:
                # Check if both reactants have aromatic rings
                aromatic_reactants = 0
                for reactant_smiles in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol:
                        aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
                        if reactant_mol.HasSubstructMatch(aromatic_pattern):
                            aromatic_reactants += 1

                # If we have at least 2 aromatic reactants and the product has an aromatic amine
                if aromatic_reactants >= 2:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol:
                        aromatic_amine_pattern = Chem.MolFromSmarts("c[NH]c")
                        if product_mol.HasSubstructMatch(aromatic_amine_pattern):
                            print(
                                f"Late-stage aromatic amine coupling detected at depth {depth}: {rsmi}"
                            )
                            late_stage_coupling = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return late_stage_coupling

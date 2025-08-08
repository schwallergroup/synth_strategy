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
    This function detects late-stage functional group modification without
    changing the carbon skeleton.
    """
    late_stage_modification = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_modification

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert SMILES to molecules
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                if product and all(r for r in reactants):
                    # Check if the carbon skeleton is preserved
                    # This is a simplified check - in reality, you would need to compare
                    # the carbon frameworks between reactants and product

                    # Check for functional group patterns
                    fg_patterns = [
                        Chem.MolFromSmarts("[C$(C=O)][NH2]"),  # Amide
                        Chem.MolFromSmarts("[C$(C=O)][O][C]"),  # Ester
                        Chem.MolFromSmarts("[OH]"),  # Alcohol
                        Chem.MolFromSmarts("[C$(C=O)]"),  # Carbonyl
                    ]

                    reactant_fgs = set()
                    product_fgs = set()

                    for i, pattern in enumerate(fg_patterns):
                        for r in reactants:
                            if r.HasSubstructMatch(pattern):
                                reactant_fgs.add(i)

                        if product.HasSubstructMatch(pattern):
                            product_fgs.add(i)

                    # If functional groups changed but we're in a late stage reaction
                    if reactant_fgs != product_fgs:
                        print("Detected late-stage functional group modification")
                        late_stage_modification = True

        # Continue traversing with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return late_stage_modification

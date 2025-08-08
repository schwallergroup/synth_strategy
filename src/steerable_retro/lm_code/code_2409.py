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
    This function detects if the synthetic route involves late-stage arylation via C-N bond formation.
    """
    late_stage_arylation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_arylation

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only check reactions at depth 0 or 1 (late stage)
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for aryl halide in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("c-[#9,#17,#35,#53]")

            # Check for new aryl-N bond in product
            aryl_n_pattern = Chem.MolFromSmarts("c-[#7;!$(N=*)]")

            reactant_has_aryl_halide = any(
                r is not None and r.HasSubstructMatch(aryl_halide_pattern) for r in reactants
            )

            if reactant_has_aryl_halide and product is not None:
                # Count aryl-N bonds in reactants and product
                reactant_aryl_n_count = sum(
                    r.GetSubstructMatches(aryl_n_pattern).__len__()
                    for r in reactants
                    if r is not None
                )
                product_aryl_n_count = product.GetSubstructMatches(aryl_n_pattern).__len__()

                if product_aryl_n_count > reactant_aryl_n_count:
                    print(f"Detected late-stage arylation via C-N bond formation at depth {depth}")
                    late_stage_arylation = True

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return late_stage_arylation

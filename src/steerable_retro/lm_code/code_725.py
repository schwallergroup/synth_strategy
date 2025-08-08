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
    This function detects if the synthetic route involves a late-stage amide coupling
    (in the first half of the synthesis depth).
    """
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    acid_chloride_pattern = Chem.MolFromSmarts("C(=O)Cl")
    amine_pattern = Chem.MolFromSmarts("[NH2]")

    max_depth = 0
    amide_coupling_depth = None

    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            find_max_depth(child, current_depth + 1)

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Check for acid chloride and amine in reactants
                    reactants_have_acid_chloride = any(
                        Chem.MolFromSmiles(r)
                        and Chem.MolFromSmiles(r).HasSubstructMatch(acid_chloride_pattern)
                        for r in reactants_smiles
                        if Chem.MolFromSmiles(r)
                    )

                    reactants_have_amine = any(
                        Chem.MolFromSmiles(r)
                        and Chem.MolFromSmiles(r).HasSubstructMatch(amine_pattern)
                        for r in reactants_smiles
                        if Chem.MolFromSmiles(r)
                    )

                    # Check for amide in product
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if (
                        reactants_have_acid_chloride
                        and reactants_have_amine
                        and product_mol
                        and product_mol.HasSubstructMatch(amide_pattern)
                    ):
                        amide_coupling_depth = depth
                        print(f"Detected amide coupling at depth {depth}")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    find_max_depth(route)
    dfs_traverse(route)

    # Check if amide coupling occurs in the first half of the synthesis
    if amide_coupling_depth is not None:
        is_late_stage = amide_coupling_depth <= (max_depth / 2)
        return is_late_stage

    return False

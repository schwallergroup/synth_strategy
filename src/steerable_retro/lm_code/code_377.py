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
    This function detects a synthetic strategy involving early amide formation
    using an acid chloride and a primary amine.
    """
    acid_chloride_pattern = Chem.MolFromSmarts("[#6](=[#8])[Cl]")
    primary_amine_pattern = Chem.MolFromSmarts("[#7;H2][#6]")
    amide_pattern = Chem.MolFromSmarts("[#6][#7][#6](=[#8])[#6]")

    found_early_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_early_amide_formation

        if (
            node["type"] == "reaction" and depth >= 2
        ):  # Only check reactions at depth 2+ (early stage)
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                if product and all(r for r in reactants):
                    # Check if reactants include acid chloride and primary amine
                    has_acid_chloride = any(
                        r.HasSubstructMatch(acid_chloride_pattern) for r in reactants
                    )
                    has_primary_amine = any(
                        r.HasSubstructMatch(primary_amine_pattern) for r in reactants
                    )

                    # Check if product has amide
                    product_has_amide = product.HasSubstructMatch(amide_pattern)

                    if has_acid_chloride and has_primary_amine and product_has_amide:
                        print(f"Found early amide formation with acid chloride: {rsmi}")
                        found_early_amide_formation = True
            except:
                print(f"Error processing reaction SMILES: {rsmi}")

        # Process children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_early_amide_formation

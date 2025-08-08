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
    Detects a sequence where an ester is hydrolyzed to a carboxylic acid,
    which is then used in an amide coupling reaction.
    """
    # Track key steps
    ester_hydrolysis_depth = None
    amide_coupling_depth = None

    # SMARTS patterns
    ester_pattern = "[#6](=[#8])-[#8]-[#6]"
    carboxylic_acid_pattern = "[#6](=[#8])-[#8;H1]"
    amide_pattern = "[#6](=[#8])-[#7]"

    def dfs_traverse(node, depth=0):
        nonlocal ester_hydrolysis_depth, amide_coupling_depth

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester hydrolysis
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts(carboxylic_acid_pattern)
                ):
                    # Check if reactants contain ester
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            Chem.MolFromSmarts(ester_pattern)
                        ):
                            ester_hydrolysis_depth = depth
                            print(f"Detected ester hydrolysis at depth {depth}")

                # Check for amide coupling
                if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts(amide_pattern)):
                    # Check if reactants contain carboxylic acid
                    has_acid = any(
                        Chem.MolFromSmiles(r).HasSubstructMatch(
                            Chem.MolFromSmarts(carboxylic_acid_pattern)
                        )
                        for r in reactants
                        if Chem.MolFromSmiles(r)
                    )
                    if has_acid:
                        amide_coupling_depth = depth
                        print(f"Detected amide coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the sequence is present (ester hydrolysis followed by amide coupling)
    has_sequence = (
        ester_hydrolysis_depth is not None
        and amide_coupling_depth is not None
        and ester_hydrolysis_depth > amide_coupling_depth
    )  # Remember: higher depth = earlier in synthesis

    if has_sequence:
        print("Detected ester hydrolysis followed by amide coupling sequence")

    return has_sequence

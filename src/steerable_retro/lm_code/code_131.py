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
    Detects a sequential functional group interconversion strategy:
    alcohol → azide → amine → sulfonamide while preserving a complex core structure.
    """
    # Initialize tracking variables
    has_alcohol_to_azide = False
    has_azide_to_amine = False
    has_amine_to_sulfonamide = False

    # SMARTS patterns for functional group detection
    alcohol_pattern = Chem.MolFromSmarts("[OH][CX4]")
    azide_pattern = Chem.MolFromSmarts("[N-]=[N+]=[N]")
    primary_amine_pattern = Chem.MolFromSmarts("[NH2][CX4]")
    sulfonamide_pattern = Chem.MolFromSmarts("[#7]-[#16](=[#8])(=[#8])-[#7]")

    # Track reactions by depth
    reactions_by_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            reactions_by_depth[depth] = node

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Traverse the route to collect reactions
    dfs_traverse(route)

    # Sort reactions by depth (from early to late in synthesis)
    sorted_depths = sorted(reactions_by_depth.keys(), reverse=True)

    # Check for the sequence of transformations
    for i, depth in enumerate(sorted_depths):
        reaction = reactions_by_depth[depth]
        rsmi = reaction["metadata"]["rsmi"]

        # Split into reactants and products
        parts = rsmi.split(">")
        reactants = parts[0].split(".")
        product = parts[-1]

        # Create RDKit molecules
        product_mol = Chem.MolFromSmiles(product)
        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

        if not product_mol or not all(reactant_mols):
            continue

        # Check for alcohol to azide transformation
        if any(
            mol.HasSubstructMatch(alcohol_pattern) for mol in reactant_mols
        ) and product_mol.HasSubstructMatch(azide_pattern):
            has_alcohol_to_azide = True
            print("Found alcohol to azide transformation at depth", depth)

        # Check for azide to amine transformation
        if any(
            mol.HasSubstructMatch(azide_pattern) for mol in reactant_mols
        ) and product_mol.HasSubstructMatch(primary_amine_pattern):
            has_azide_to_amine = True
            print("Found azide to amine transformation at depth", depth)

        # Check for amine to sulfonamide transformation
        if any(
            mol.HasSubstructMatch(primary_amine_pattern) for mol in reactant_mols
        ) and product_mol.HasSubstructMatch(sulfonamide_pattern):
            has_amine_to_sulfonamide = True
            print("Found amine to sulfonamide transformation at depth", depth)

    # Check if the complete sequence was found
    sequence_found = has_alcohol_to_azide and has_azide_to_amine and has_amine_to_sulfonamide

    if sequence_found:
        print("Complete alcohol → azide → amine → sulfonamide sequence detected")
    else:
        print(
            "Incomplete sequence: alcohol→azide:",
            has_alcohol_to_azide,
            "azide→amine:",
            has_azide_to_amine,
            "amine→sulfonamide:",
            has_amine_to_sulfonamide,
        )

    return sequence_found

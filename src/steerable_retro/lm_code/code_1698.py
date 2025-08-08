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
    This function detects a late-stage Suzuki coupling strategy where a boronic acid derivative
    is prepared in early stages and used in a late-stage coupling to join complex fragments.
    """
    # Track if we found the key reactions
    found_suzuki_coupling = False
    found_borylation = False
    suzuki_depth = None
    borylation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal found_suzuki_coupling, found_borylation, suzuki_depth, borylation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and all(reactants):
                # Check for Suzuki coupling
                boronic_pattern = Chem.MolFromSmarts("[B][O]")
                halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")

                has_boronic = any(mol.HasSubstructMatch(boronic_pattern) for mol in reactants)
                has_halide = any(mol.HasSubstructMatch(halide_pattern) for mol in reactants)

                if has_boronic and has_halide and depth <= 1:  # Late-stage (depth 0 or 1)
                    found_suzuki_coupling = True
                    suzuki_depth = depth
                    print(f"Found Suzuki coupling at depth {depth}")

                # Check for borylation
                boronic_ester_pattern = Chem.MolFromSmarts("[c][B]1O[C](C)(C)[C](C)(C)O1")
                if (
                    any(mol.HasSubstructMatch(boronic_ester_pattern) for mol in reactants)
                    and depth >= 2
                ):  # Early-stage
                    found_borylation = True
                    borylation_depth = depth
                    print(f"Found borylation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found the pattern: borylation in early stages, Suzuki in late stages
    if found_suzuki_coupling and found_borylation:
        if (
            suzuki_depth is not None
            and borylation_depth is not None
            and suzuki_depth < borylation_depth
        ):
            print("Detected late-stage Suzuki coupling strategy")
            return True

    return False

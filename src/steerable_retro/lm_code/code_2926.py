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
    This function detects biaryl formation via Suzuki coupling.
    """
    boronic_acid_pattern = Chem.MolFromSmarts("[#6]B(O)(O)")
    halide_patterns = [
        Chem.MolFromSmarts("[#6][Br]"),
        Chem.MolFromSmarts("[#6][I]"),
        Chem.MolFromSmarts("[#6][Cl]"),
    ]
    suzuki_detected = False

    def dfs_traverse(node):
        nonlocal suzuki_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if reactants_mol and product_mol:
                        # Check for boronic acid in reactants
                        has_boronic_acid = reactants_mol.HasSubstructMatch(boronic_acid_pattern)

                        # Check for halide in reactants
                        has_halide = any(
                            reactants_mol.HasSubstructMatch(pattern) for pattern in halide_patterns
                        )

                        # Count aromatic rings in reactants and product
                        reactants_rings = len(Chem.GetSSSR(reactants_mol))
                        product_rings = len(Chem.GetSSSR(product_mol))

                        # If boronic acid and halide are present in reactants, and the number of rings is preserved or increased
                        # in the product, it's likely a Suzuki coupling
                        if has_boronic_acid and has_halide and product_rings >= reactants_rings:
                            suzuki_detected = True
                            print(
                                f"Biaryl formation via Suzuki coupling detected in reaction: {rsmi}"
                            )
                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_detected

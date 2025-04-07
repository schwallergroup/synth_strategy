#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects if the synthesis route involves late-stage functionalization,
    specifically chlorination and/or alkoxy group installation in the last steps.
    """
    # Track late-stage functionalizations
    late_functionalizations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product_mol and any(reactant_mols):
                # Find the main reactant (usually the most complex one)
                main_reactant = max(
                    reactant_mols, key=lambda m: m.GetNumAtoms() if m else 0
                )

                # Check for chlorination
                chlorine_pattern = Chem.MolFromSmarts("[#6]-[Cl]")
                reactant_has_chlorine = (
                    main_reactant.HasSubstructMatch(chlorine_pattern)
                    if main_reactant
                    else False
                )
                product_has_chlorine = product_mol.HasSubstructMatch(chlorine_pattern)

                if product_has_chlorine and not reactant_has_chlorine:
                    late_functionalizations.append(("chlorination", depth))
                    print(f"Chlorination at depth {depth}")

                # Check for alkoxy group installation
                alkoxy_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]")

                # Count alkoxy groups in reactant and product
                reactant_alkoxy_count = (
                    len(main_reactant.GetSubstructMatches(alkoxy_pattern))
                    if main_reactant
                    else 0
                )
                product_alkoxy_count = len(
                    product_mol.GetSubstructMatches(alkoxy_pattern)
                )

                if product_alkoxy_count > reactant_alkoxy_count:
                    late_functionalizations.append(("alkoxy_installation", depth))
                    print(f"Alkoxy group installation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Late-stage functionalization is defined as occurring at depth <= 1
    return any(depth <= 1 for _, depth in late_functionalizations)

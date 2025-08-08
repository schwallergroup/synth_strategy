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
    This function detects a synthetic strategy involving Grignard addition to an
    α,β-unsaturated system.
    """
    has_grignard_addition = False

    def dfs_traverse(node, depth=0):
        nonlocal has_grignard_addition

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Look for patterns suggesting Grignard addition
            for reactant in reactants:
                if "MgBr" in reactant or "BrMg" in reactant:
                    # Check if other reactant has α,β-unsaturated system
                    for other_reactant in reactants:
                        if other_reactant != reactant:
                            other_mol = Chem.MolFromSmiles(other_reactant)
                            if other_mol:
                                # α,β-unsaturated system pattern
                                unsaturated_pattern = Chem.MolFromSmarts("[#6]=[#6][#6](=[#8,#7])")
                                if other_mol.HasSubstructMatch(unsaturated_pattern):
                                    has_grignard_addition = True
                                    print(f"Found Grignard addition at depth {depth}")

            # Alternative detection method using reaction pattern
            if len(reactants) >= 2:
                # Check if one reactant is aromatic and product has new C-C bond to aromatic
                aromatic_pattern = Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1")
                unsaturated_pattern = Chem.MolFromSmarts("[#6]=[#6][#6](=[#8,#7])")

                product_mol = Chem.MolFromSmiles(product)

                has_aromatic_reactant = False
                has_unsaturated_reactant = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(aromatic_pattern):
                            has_aromatic_reactant = True
                        if reactant_mol.HasSubstructMatch(unsaturated_pattern):
                            has_unsaturated_reactant = True

                if (
                    has_aromatic_reactant
                    and has_unsaturated_reactant
                    and product_mol
                    and product_mol.HasSubstructMatch(aromatic_pattern)
                ):
                    has_grignard_addition = True
                    print(f"Found likely Grignard addition at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Grignard addition strategy detected: {has_grignard_addition}")

    return has_grignard_addition

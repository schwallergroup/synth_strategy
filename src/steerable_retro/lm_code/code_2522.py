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
    This function detects a strategy involving multiple Friedel-Crafts acylation
    reactions on aromatic rings.
    """
    friedel_crafts_reactions = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Friedel-Crafts acylation pattern
                try:
                    product_mol = Chem.MolFromSmiles(product)

                    # Look for acyl chloride in reactants
                    acyl_chloride_present = False
                    aromatic_substrate_present = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            acyl_chloride_pattern = Chem.MolFromSmarts("[C](=[O])[Cl]")
                            if reactant_mol.HasSubstructMatch(acyl_chloride_pattern):
                                acyl_chloride_present = True

                            aromatic_pattern = Chem.MolFromSmarts("c")
                            if reactant_mol.HasSubstructMatch(aromatic_pattern):
                                aromatic_substrate_present = True

                    # Check if product has new ketone attached to aromatic ring
                    if acyl_chloride_present and aromatic_substrate_present:
                        aromatic_ketone_pattern = Chem.MolFromSmarts("c[C](=[O])[#6]")
                        if product_mol.HasSubstructMatch(aromatic_ketone_pattern):
                            depth = node.get("depth", "unknown")
                            friedel_crafts_reactions.append(depth)
                            print(f"Found Friedel-Crafts acylation at depth: {depth}")

                except Exception as e:
                    print(f"Error in Friedel-Crafts detection: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least 2 Friedel-Crafts reactions
    has_sequential_fc = len(friedel_crafts_reactions) >= 2
    print(f"Friedel-Crafts reactions found at depths: {friedel_crafts_reactions}")
    return has_sequential_fc

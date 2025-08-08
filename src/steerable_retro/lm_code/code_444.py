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
    This function detects a synthetic strategy where an aromatic core structure
    (particularly dichlorobenzene) is preserved throughout the synthesis while
    side chains are modified.
    """
    # Track molecules in the main synthetic pathway
    main_pathway_molecules = []

    def dfs_traverse(node, is_main_product=True):
        if node["type"] == "mol":
            # Only track molecules in the main synthetic pathway
            if is_main_product and not node.get("in_stock", False):
                main_pathway_molecules.append(node["smiles"])
                print(f"Adding to main pathway: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # For reaction nodes, only the product is in the main pathway
            # Children are reactants in retrosynthetic direction
            for i, child in enumerate(node.get("children", [])):
                # First child is typically the main reactant in retrosynthesis
                is_main_reactant = i == 0
                dfs_traverse(child, is_main_reactant)
            return  # Skip the regular children traversal

        # Process children for non-reaction nodes
        for child in node.get("children", []):
            dfs_traverse(child, is_main_product)

    # Traverse the route to collect main pathway molecules
    dfs_traverse(route)

    # Check if all main pathway molecules have the dichlorobenzene core
    core_preserved = True
    for smiles in main_pathway_molecules:
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.GetNumAtoms() > 6:  # Skip small molecules
            # Check for dichlorobenzene pattern
            dichlorobenzene_pattern = Chem.MolFromSmarts("c1c(Cl)cc(Cl)cc1")
            if mol.HasSubstructMatch(dichlorobenzene_pattern):
                print(f"Core preserved in molecule: {smiles}")
            else:
                print(f"Core not preserved in molecule: {smiles}")
                core_preserved = False
                break

    return core_preserved

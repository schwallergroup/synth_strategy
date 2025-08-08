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
    Detects if the synthesis route uses cross-coupling reactions (like Suzuki coupling)
    for scaffold construction, particularly for C-C bond formation between aromatic systems.
    """
    cross_coupling_detected = False

    def dfs_traverse(node):
        nonlocal cross_coupling_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronic acid pattern in reactants
                boronic_acid_pattern = Chem.MolFromSmarts("[#6]-[B]([O])([O])")

                # Check for halogen pattern in reactants
                halogen_pattern = Chem.MolFromSmarts("[#6]-[#53,#35,#17]")  # C-I, C-Br, C-Cl

                # Convert SMILES to molecules
                try:
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                    # Check if one reactant has boronic acid and another has halogen
                    has_boronic_acid = any(
                        mol and mol.HasSubstructMatch(boronic_acid_pattern) for mol in reactant_mols
                    )
                    has_halogen = any(
                        mol and mol.HasSubstructMatch(halogen_pattern) for mol in reactant_mols
                    )

                    if has_boronic_acid and has_halogen:
                        print(
                            "Cross-coupling reaction detected: Boronic acid and halogen reactants found"
                        )
                        cross_coupling_detected = True
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return cross_coupling_detected

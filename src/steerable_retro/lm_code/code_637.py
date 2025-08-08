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
    This function detects a synthetic strategy involving coupling of two aromatic fragments
    (thioanisole and fluorophenyl-containing fragments).
    """
    thioanisole_fragment = False
    fluorophenyl_fragment = False
    coupling_reaction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal thioanisole_fragment, fluorophenyl_fragment, coupling_reaction_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert SMILES to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product_mol and len(reactant_mols) > 1:
                # Check for thioanisole and fluorophenyl fragments
                thioanisole_pattern = Chem.MolFromSmarts("c-S-[CH3]")
                fluorophenyl_pattern = Chem.MolFromSmarts("c-F")

                thioanisole_in_reactants = any(
                    r.HasSubstructMatch(thioanisole_pattern) for r in reactant_mols
                )
                fluorophenyl_in_reactants = any(
                    r.HasSubstructMatch(fluorophenyl_pattern) for r in reactant_mols
                )

                if (
                    thioanisole_in_reactants
                    and fluorophenyl_in_reactants
                    and product_mol.HasSubstructMatch(thioanisole_pattern)
                    and product_mol.HasSubstructMatch(fluorophenyl_pattern)
                ):
                    print(f"Aromatic fragment coupling detected at depth {depth}")
                    thioanisole_fragment = True
                    fluorophenyl_fragment = True
                    coupling_reaction_found = True

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    strategy_present = thioanisole_fragment and fluorophenyl_fragment and coupling_reaction_found
    print(f"Aromatic fragment coupling strategy detected: {strategy_present}")
    return strategy_present

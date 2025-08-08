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
    Detects a strategy involving heterocyclic ring formation (specifically benzoxazole)
    as a key step in the synthesis.
    """
    heterocycle_formed = False

    def dfs_traverse(node):
        nonlocal heterocycle_formed

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]

                # Check for benzoxazole formation
                benzoxazole_pattern = Chem.MolFromSmarts("c1nc(C)oc1")

                if product_mol and product_mol.HasSubstructMatch(benzoxazole_pattern):
                    reactant_has_benzoxazole = any(
                        r and r.HasSubstructMatch(benzoxazole_pattern) for r in reactants_mols if r
                    )

                    if not reactant_has_benzoxazole:
                        # Check if reactants contain phenol and acetamide patterns
                        phenol_pattern = Chem.MolFromSmarts("c[OH]")
                        acetamide_pattern = Chem.MolFromSmarts("[NH]C(=O)C")

                        reactant_has_phenol = any(
                            r and r.HasSubstructMatch(phenol_pattern) for r in reactants_mols if r
                        )

                        reactant_has_acetamide = any(
                            r and r.HasSubstructMatch(acetamide_pattern)
                            for r in reactants_mols
                            if r
                        )

                        if reactant_has_phenol or reactant_has_acetamide:
                            print(
                                f"Benzoxazole formation from phenol/acetamide detected at depth: {node.get('metadata', {}).get('depth', 'unknown')}"
                            )
                            heterocycle_formed = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return heterocycle_formed

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
    This function detects if the synthetic route involves N-alkylation of a heterocycle
    in the middle stages of the synthesis.
    """
    n_alkylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            try:
                # Look for halogenated alkyl groups in reactants
                alkyl_halide_pattern = Chem.MolFromSmarts("[#6]-[#9,#17,#35,#53]")

                # Look for heterocyclic nitrogen in reactants
                n_heterocycle_pattern = Chem.MolFromSmarts("[#7;R]")

                # Check if product has N-C bond where N was previously in heterocycle
                n_alkyl_pattern = Chem.MolFromSmarts("[#7;R]-[#6;!R]")

                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if all(reactants) and product_mol:
                    has_alkyl_halide = any(
                        mol.HasSubstructMatch(alkyl_halide_pattern) for mol in reactants
                    )
                    has_n_heterocycle = any(
                        mol.HasSubstructMatch(n_heterocycle_pattern) for mol in reactants
                    )

                    if (
                        has_alkyl_halide
                        and has_n_heterocycle
                        and product_mol.HasSubstructMatch(n_alkyl_pattern)
                    ):
                        print(f"Found N-alkylation of heterocycle at depth {depth}")
                        n_alkylation_found = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return n_alkylation_found

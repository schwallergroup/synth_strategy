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
    This function detects synthesis routes that involve SNAr coupling (aromatic nucleophilic substitution)
    where a halogen is replaced by a nitrogen nucleophile.
    """
    has_snar_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_snar_coupling

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for SNAr coupling
            # Look for aromatic C-N bonds in product
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol:
                # Find aromatic C-N bonds in product
                aromatic_c_n_pattern = Chem.MolFromSmarts("c[N]")
                if aromatic_c_n_pattern and product_mol.HasSubstructMatch(aromatic_c_n_pattern):
                    # Check if any reactant has a halogen on aromatic carbon
                    halogen_aromatic_pattern = Chem.MolFromSmarts("c[F,Cl,Br,I]")
                    for r in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r)
                        if (
                            r_mol
                            and halogen_aromatic_pattern
                            and r_mol.HasSubstructMatch(halogen_aromatic_pattern)
                        ):
                            # Check if another reactant has an amine
                            amine_pattern = Chem.MolFromSmarts("[NH2]")
                            for r2 in reactants_smiles:
                                if r2 != r:
                                    r2_mol = Chem.MolFromSmiles(r2)
                                    if (
                                        r2_mol
                                        and amine_pattern
                                        and r2_mol.HasSubstructMatch(amine_pattern)
                                    ):
                                        has_snar_coupling = True
                                        print(f"Found SNAr coupling at depth {depth}")
                                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_snar_coupling

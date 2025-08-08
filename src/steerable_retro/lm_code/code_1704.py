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
    Detects a synthetic strategy where SNAr occurs on an activated aromatic ring
    with electron-withdrawing groups (nitro, fluoro) ortho or para to the leaving group.
    """
    # Track if we found an activated SNAr
    found_activated_snar = False

    def dfs_traverse(node, depth=0):
        nonlocal found_activated_snar

        if node["type"] == "reaction":
            # Extract reactants and product
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for potential SNAr reaction
                # Look for halogen in reactants
                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    # Check for aryl halide with activating groups
                    # Ortho activation patterns
                    ortho_nitro_halide = Chem.MolFromSmarts("c([N+](=O)[O-])c[F,Cl,Br,I]")
                    ortho_fluoro_halide = Chem.MolFromSmarts("c(F)c[F,Cl,Br,I]")

                    # Para activation patterns
                    para_nitro_halide = Chem.MolFromSmarts("c([N+](=O)[O-])ccc[F,Cl,Br,I]")
                    para_fluoro_halide = Chem.MolFromSmarts("c(F)ccc[F,Cl,Br,I]")

                    if (
                        reactant_mol.HasSubstructMatch(ortho_nitro_halide)
                        or reactant_mol.HasSubstructMatch(ortho_fluoro_halide)
                        or reactant_mol.HasSubstructMatch(para_nitro_halide)
                        or reactant_mol.HasSubstructMatch(para_fluoro_halide)
                    ):

                        # Check if product has a new C-N bond where the halogen was
                        product_mol = Chem.MolFromSmiles(product_smiles)
                        if product_mol:
                            # Look for aromatic C-N bond that wasn't in the reactant
                            c_n_bond_pattern = Chem.MolFromSmarts("cN")
                            if product_mol.HasSubstructMatch(c_n_bond_pattern):
                                print(f"Found activated SNAr at depth {depth}")
                                found_activated_snar = True

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_activated_snar

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
    This function detects a strategy involving sequential heteroatom alkylations on a beta-lactam core,
    specifically an early S-alkylation followed by a late N-alkylation.
    """
    # Track if we found the required reactions
    found_s_alkylation = False
    found_n_alkylation = False
    s_alkylation_depth = -1
    n_alkylation_depth = -1

    # Beta-lactam pattern
    beta_lactam_pattern = Chem.MolFromSmarts("[#6]1[#7][#6](=[#8])[#6]1")

    def dfs_traverse(node, depth=0):
        nonlocal found_s_alkylation, found_n_alkylation, s_alkylation_depth, n_alkylation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Check if all molecules were parsed correctly
            if None in reactant_mols or product_mol is None:
                print("Warning: Could not parse some molecules in reaction")
                return

            # Check if product contains beta-lactam
            if product_mol.HasSubstructMatch(beta_lactam_pattern):
                # Check for S-alkylation (C-S bond formation)
                if len(reactant_mols) == 2:  # Bimolecular reaction
                    # Look for SH in one reactant and CH2-X in another
                    sh_pattern = Chem.MolFromSmarts("[#6]-[#16;H1]")
                    benzyl_halide_pattern = Chem.MolFromSmarts("[#6]-[#6;H2]-[#9,#17,#35,#53]")

                    r1_has_sh = any(r.HasSubstructMatch(sh_pattern) for r in reactant_mols)
                    r2_has_benzyl_halide = any(
                        r.HasSubstructMatch(benzyl_halide_pattern) for r in reactant_mols
                    )

                    # Check if product has C-S-C pattern (thioether)
                    thioether_pattern = Chem.MolFromSmarts("[#6]-[#16]-[#6]")
                    product_has_thioether = product_mol.HasSubstructMatch(thioether_pattern)

                    if r1_has_sh and r2_has_benzyl_halide and product_has_thioether:
                        found_s_alkylation = True
                        s_alkylation_depth = depth
                        print(f"Found S-alkylation at depth {depth}")

                # Check for N-alkylation of beta-lactam
                lactam_n_pattern = Chem.MolFromSmarts(
                    "[#6]1[#7;H1][#6](=[#8])[#6]1"
                )  # NH in beta-lactam
                alkylated_lactam_pattern = Chem.MolFromSmarts(
                    "[#6]1[#7;!H1][#6](=[#8])[#6]1"
                )  # N-substituted beta-lactam

                r_has_nh_lactam = any(r.HasSubstructMatch(lactam_n_pattern) for r in reactant_mols)
                product_has_n_alkylated_lactam = product_mol.HasSubstructMatch(
                    alkylated_lactam_pattern
                )

                if r_has_nh_lactam and product_has_n_alkylated_lactam:
                    found_n_alkylation = True
                    n_alkylation_depth = depth
                    print(f"Found N-alkylation at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found both reactions in the correct order (S-alkylation early, N-alkylation late)
    if found_s_alkylation and found_n_alkylation and s_alkylation_depth > n_alkylation_depth:
        print("Found beta-lactam sequential heteroatom alkylation strategy")
        return True

    return False

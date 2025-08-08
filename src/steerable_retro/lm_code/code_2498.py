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
    Detects a synthetic strategy where an early Suzuki coupling creates a biaryl scaffold,
    followed by nitro reduction, amine functionalization, and final amide formation.
    """
    # Initialize tracking variables
    has_biaryl_formation = False
    has_nitro_reduction = False
    has_amine_protection = False
    has_amine_alkylation = False
    has_amide_formation = False
    biaryl_depth = -1
    nitro_reduction_depth = -1
    amide_formation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal has_biaryl_formation, has_nitro_reduction, has_amine_protection
        nonlocal has_amine_alkylation, has_amide_formation
        nonlocal biaryl_depth, nitro_reduction_depth, amide_formation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product_mol and len(reactant_mols) >= 1:
                # Check for Suzuki coupling (biaryl formation)
                if len(reactant_mols) >= 2:
                    # Look for boronic acid/ester and aryl halide patterns
                    boronic_acid_pattern = Chem.MolFromSmarts("[c]-[B]")
                    aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Cl,Br,I,F]")

                    has_boronic_acid = any(
                        mol.HasSubstructMatch(boronic_acid_pattern) for mol in reactant_mols if mol
                    )
                    has_aryl_halide = any(
                        mol.HasSubstructMatch(aryl_halide_pattern) for mol in reactant_mols if mol
                    )

                    # Check if product has new biaryl bond
                    biaryl_pattern = Chem.MolFromSmarts("c-c")

                    if (
                        has_boronic_acid
                        and has_aryl_halide
                        and product_mol.HasSubstructMatch(biaryl_pattern)
                    ):
                        has_biaryl_formation = True
                        biaryl_depth = depth
                        print(f"Detected biaryl formation at depth {depth}")

                # Check for nitro reduction
                nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])-[O-]")
                amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")

                has_nitro_reactant = any(
                    mol.HasSubstructMatch(nitro_pattern) for mol in reactant_mols if mol
                )
                has_amine_product = product_mol and product_mol.HasSubstructMatch(amine_pattern)

                if has_nitro_reactant and has_amine_product:
                    has_nitro_reduction = True
                    nitro_reduction_depth = depth
                    print(f"Detected nitro reduction at depth {depth}")

                # Check for Boc protection
                boc_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[O])-[#7]")
                amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")
                protected_amine_pattern = Chem.MolFromSmarts("[#6]-[#7]-[#6](=[O])-[#8]-[#6]")

                has_amine_reactant = any(
                    mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols if mol
                )
                has_protected_amine_product = product_mol and product_mol.HasSubstructMatch(
                    protected_amine_pattern
                )

                if has_amine_reactant and has_protected_amine_product:
                    has_amine_protection = True
                    print(f"Detected amine protection at depth {depth}")

                # Check for amine alkylation (focusing on methylation)
                amine_pattern = Chem.MolFromSmarts("[#6]-[NH]-[#6]")
                has_methylamine_product = product_mol and product_mol.HasSubstructMatch(
                    amine_pattern
                )

                if has_methylamine_product:
                    has_amine_alkylation = True
                    print(f"Detected amine alkylation at depth {depth}")

                # Check for amide formation
                carboxylic_acid_pattern = Chem.MolFromSmarts("[#6]-[#6](=[O])-[#8]")
                amine_pattern = Chem.MolFromSmarts("[#6]-[#7](-[#6])-[#1]")
                amide_pattern = Chem.MolFromSmarts("[#6]-[#6](=[O])-[#7](-[#6])-[#6]")

                has_acid_reactant = any(
                    mol.HasSubstructMatch(carboxylic_acid_pattern) for mol in reactant_mols if mol
                )
                has_amine_reactant = any(
                    mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols if mol
                )
                has_amide_product = product_mol and product_mol.HasSubstructMatch(amide_pattern)

                if (has_acid_reactant or has_amine_reactant) and has_amide_product:
                    has_amide_formation = True
                    amide_formation_depth = depth
                    print(f"Detected amide formation at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        has_biaryl_formation
        and has_nitro_reduction
        and (has_amine_protection or has_amine_alkylation)
        and has_amide_formation
        and biaryl_depth > nitro_reduction_depth > amide_formation_depth
    )

    if strategy_present:
        print("Detected strategy: Biaryl formation followed by amine functionalization")

    return strategy_present

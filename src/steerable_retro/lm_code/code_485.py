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
    Detects if the synthesis involves sequential C-N bond formations
    (amide formation, acylation, nucleophilic substitution).
    """
    # Track C-N bond formations
    cn_bond_formations = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for amide formation
            amide_pattern = Chem.MolFromSmarts("[NH][C](=[O])")
            amine_pattern = Chem.MolFromSmarts("[NH2]")
            acyl_pattern = Chem.MolFromSmarts("[C](=[O])[Cl]")

            # Check for nucleophilic substitution
            alkyl_halide_pattern = Chem.MolFromSmarts("[CH2][Cl]")
            amine_nucleophile_pattern = Chem.MolFromSmarts("[N]1[CH2][CH2][N][CH2][CH2]1")

            has_amine_reactant = False
            has_acyl_reactant = False
            has_alkyl_halide = False
            has_amine_nucleophile = False

            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if not mol:
                        continue

                    if mol.HasSubstructMatch(amine_pattern):
                        has_amine_reactant = True
                    if mol.HasSubstructMatch(acyl_pattern):
                        has_acyl_reactant = True
                    if mol.HasSubstructMatch(alkyl_halide_pattern):
                        has_alkyl_halide = True
                    if mol.HasSubstructMatch(amine_nucleophile_pattern):
                        has_amine_nucleophile = True
                except:
                    continue

            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                    if has_amine_reactant and has_acyl_reactant:
                        cn_bond_formations.append("amide_formation")
                        print("Detected amide formation")

                # Check for nucleophilic substitution product
                if has_alkyl_halide and has_amine_nucleophile:
                    cn_bond_formations.append("nucleophilic_substitution")
                    print("Detected nucleophilic substitution")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if we have at least 2 different types of C-N bond formations
    return len(set(cn_bond_formations)) >= 2

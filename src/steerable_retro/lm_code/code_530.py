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
    Detects a sequential alkylation strategy on a piperidone scaffold:
    1. C-alkylation with halogenated aryl group
    2. Decarboxylation step
    3. N-alkylation with benzyl group
    """
    # Track if we found each step
    found_c_alkylation = False
    found_decarboxylation = False
    found_n_alkylation = False

    # Track reaction depths for ordering
    c_alkylation_depth = -1
    decarboxylation_depth = -1
    n_alkylation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_c_alkylation, found_decarboxylation, found_n_alkylation
        nonlocal c_alkylation_depth, decarboxylation_depth, n_alkylation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
            product = Chem.MolFromSmiles(product_smiles)

            if product is None or any(r is None for r in reactants):
                print(f"Warning: Could not parse SMILES at depth {depth}")
                return

            # Check for piperidone scaffold in product
            piperidone_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#7][#6](=[O])1")
            if product.HasSubstructMatch(piperidone_pattern):

                # Check for C-alkylation with halogenated aryl group
                if len(reactants) >= 2:
                    # Look for dichlorobenzyl pattern in one of the reactants
                    dichlorobenzyl_pattern = Chem.MolFromSmarts(
                        "[#6]-[c]1[c]([Cl])[c][c][c][c]1[Cl]"
                    )
                    alkyl_halide_pattern = Chem.MolFromSmarts("[C]-[Br,Cl,I,F]")

                    for r in reactants:
                        if r.HasSubstructMatch(dichlorobenzyl_pattern) and r.HasSubstructMatch(
                            alkyl_halide_pattern
                        ):
                            # Check if this is forming a C-C bond to the piperidone
                            if "[C:6]1([CH2:7]" in rsmi or "[CH:6]1[CH2:7]" in rsmi:
                                found_c_alkylation = True
                                c_alkylation_depth = depth
                                print(f"Found C-alkylation with halogenated aryl at depth {depth}")

                # Check for N-alkylation with benzyl group
                if len(reactants) >= 2:
                    # Look for benzyl pattern in one of the reactants
                    benzyl_pattern = Chem.MolFromSmarts("[#6]-[c]1[c][c][c][c][c]1")
                    alkyl_halide_pattern = Chem.MolFromSmarts("[C]-[Br,Cl,I,F]")

                    for r in reactants:
                        if r.HasSubstructMatch(benzyl_pattern) and r.HasSubstructMatch(
                            alkyl_halide_pattern
                        ):
                            # Check if this is forming a C-N bond to the piperidone
                            if "[NH:16]1" in reactants_smiles and "[N:16]1" in product_smiles:
                                found_n_alkylation = True
                                n_alkylation_depth = depth
                                print(f"Found N-alkylation with benzyl at depth {depth}")

                # Check for decarboxylation
                ester_pattern = Chem.MolFromSmarts("[C]-[C](=[O])-[O]-[C]")
                if any(
                    r.HasSubstructMatch(ester_pattern) for r in reactants
                ) and not product.HasSubstructMatch(ester_pattern):
                    found_decarboxylation = True
                    decarboxylation_depth = depth
                    print(f"Found decarboxylation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found all steps in the correct order
    correct_order = (
        found_c_alkylation
        and found_decarboxylation
        and found_n_alkylation
        and c_alkylation_depth > decarboxylation_depth > n_alkylation_depth
    )

    if correct_order:
        print("Found sequential alkylation strategy on piperidone scaffold")
    else:
        print("Did not find complete sequential alkylation strategy")
        print(f"C-alkylation: {found_c_alkylation} at depth {c_alkylation_depth}")
        print(f"Decarboxylation: {found_decarboxylation} at depth {decarboxylation_depth}")
        print(f"N-alkylation: {found_n_alkylation} at depth {n_alkylation_depth}")

    return correct_order

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
    This function detects a protect-couple-deprotect sequence for carboxylic acids.
    Looks for carboxylic acid → ester → coupling → carboxylic acid sequence.
    """
    # Track steps in the sequence
    esterification_step = None
    coupling_step = None
    deprotection_step = None
    step_depths = {}

    def dfs_traverse(node, depth=0):
        nonlocal esterification_step, coupling_step, deprotection_step

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for esterification (carboxylic acid → ester)
            acid_pattern = Chem.MolFromSmarts("[OH][C]=O")
            ester_pattern = Chem.MolFromSmarts("[#6][O][C]=O")

            # Check for coupling reaction (look for C-C bond formation)
            coupling_pattern = Chem.MolFromSmarts("[c]!@[c]")

            # Check for deprotection (ester → carboxylic acid)

            try:
                # Check for esterification
                has_acid_reactant = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(acid_pattern):
                        has_acid_reactant = True
                        break

                product_mol = Chem.MolFromSmiles(product)
                has_ester_product = product_mol and product_mol.HasSubstructMatch(ester_pattern)

                if has_acid_reactant and has_ester_product:
                    esterification_step = depth
                    step_depths[depth] = "esterification"
                    print(f"Esterification detected at depth {depth}")

                # Check for coupling
                if product_mol and product_mol.HasSubstructMatch(coupling_pattern):
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and not mol.HasSubstructMatch(coupling_pattern):
                            coupling_step = depth
                            step_depths[depth] = "coupling"
                            print(f"Coupling reaction detected at depth {depth}")
                            break

                # Check for deprotection
                has_ester_reactant = False
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(ester_pattern):
                        has_ester_reactant = True
                        break

                has_acid_product = product_mol and product_mol.HasSubstructMatch(acid_pattern)

                if has_ester_reactant and has_acid_product:
                    deprotection_step = depth
                    step_depths[depth] = "deprotection"
                    print(f"Deprotection detected at depth {depth}")

            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have all three steps in the correct order
    if (
        esterification_step is not None
        and coupling_step is not None
        and deprotection_step is not None
    ):
        # In retrosynthetic direction, deprotection should be earlier (lower depth) than coupling,
        # and coupling should be earlier than esterification
        if deprotection_step < coupling_step < esterification_step:
            print("Protect-couple-deprotect sequence detected for carboxylic acid")
            return True

    return False

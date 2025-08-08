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
    This function detects a protection-deprotection sequence in the synthetic route,
    specifically looking for ketone protection via ketal formation.
    """
    ketone_pattern = Chem.MolFromSmarts("[#6][C](=[O])[#6]")
    ketal_pattern = Chem.MolFromSmarts("[#6]1[#8][#6][#6][#8]1")
    protection_step = False
    deprotection_step = False

    def dfs_traverse(node):
        nonlocal protection_step, deprotection_step

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                # Check for protection: ketone + diol -> ketal
                if product_mol and product_mol.HasSubstructMatch(ketal_pattern):
                    has_ketone = any(
                        r_mol and r_mol.HasSubstructMatch(ketone_pattern)
                        for r_mol in reactants_mols
                        if r_mol
                    )
                    if has_ketone:
                        protection_step = True
                        print("Detected ketone protection step")

                # Check for deprotection: ketal -> ketone + diol
                has_ketal = any(
                    r_mol and r_mol.HasSubstructMatch(ketal_pattern)
                    for r_mol in reactants_mols
                    if r_mol
                )
                if has_ketal and product_mol and product_mol.HasSubstructMatch(ketone_pattern):
                    deprotection_step = True
                    print("Detected ketal deprotection step")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both protection and deprotection steps are detected
    return protection_step and deprotection_step

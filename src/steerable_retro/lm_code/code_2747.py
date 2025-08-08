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
    This function detects a sequence of aromatic functionalization steps,
    specifically looking for O-alkylation and halogenation on aromatic rings.
    """
    # Track functionalization steps
    o_alkylation_detected = False
    halogenation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal o_alkylation_detected, halogenation_detected

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if all(m is not None for m in reactant_mols) and product_mol is not None:
                    # Check for O-alkylation (phenol to methoxy)
                    phenol_pattern = Chem.MolFromSmarts("[#8H][#6]:[#6]")
                    methoxy_pattern = Chem.MolFromSmarts("[#8][#6][#6]:[#6]")

                    if any(
                        mol.HasSubstructMatch(phenol_pattern) for mol in reactant_mols
                    ) and product_mol.HasSubstructMatch(methoxy_pattern):
                        o_alkylation_detected = True
                        print(f"O-alkylation detected at depth {depth}")

                    # Check for aromatic halogenation
                    reactant_bromo_count = sum(
                        len(mol.GetSubstructMatches(Chem.MolFromSmarts("c[Br]")))
                        for mol in reactant_mols
                    )
                    product_bromo_count = len(
                        product_mol.GetSubstructMatches(Chem.MolFromSmarts("c[Br]"))
                    )

                    if product_bromo_count > reactant_bromo_count:
                        halogenation_detected = True
                        print(f"Aromatic halogenation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if both functionalization types were detected
    if o_alkylation_detected and halogenation_detected:
        print("Sequential aromatic functionalization strategy detected")
        return True
    return False

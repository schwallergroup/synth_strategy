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
    Detects if the synthesis route involves a sequence of aromatic functionalizations.
    Looks for multiple modifications to an aromatic ring.
    """
    aromatic_modifications = 0

    def dfs_traverse(node):
        nonlocal aromatic_modifications

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Patterns for common aromatic modifications
            nitration_pattern = Chem.MolFromSmarts("[c][N+](=[O])[O-]")
            halogenation_pattern = Chem.MolFromSmarts("[c][F,Cl,Br,I]")
            alkoxy_pattern = Chem.MolFromSmarts("[c][O][C]")
            borylation_pattern = Chem.MolFromSmarts("[c][B]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            # Check if the reaction modifies an aromatic ring
            if product_mol is not None:
                for pattern in [
                    nitration_pattern,
                    halogenation_pattern,
                    alkoxy_pattern,
                    borylation_pattern,
                ]:
                    product_matches = (
                        product_mol.GetSubstructMatches(pattern)
                        if product_mol.HasSubstructMatch(pattern)
                        else []
                    )
                    reactant_matches = sum(
                        [
                            (
                                len(mol.GetSubstructMatches(pattern))
                                if mol is not None and mol.HasSubstructMatch(pattern)
                                else 0
                            )
                            for mol in reactant_mols
                        ]
                    )

                    if len(product_matches) != reactant_matches:
                        print(
                            f"Found aromatic modification at depth {node.get('metadata', {}).get('ID', '')}"
                        )
                        aromatic_modifications += 1
                        break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return aromatic_modifications >= 3  # At least 3 modifications to count as a sequence

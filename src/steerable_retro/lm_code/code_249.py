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
    Detects if the synthetic route contains a fragment coupling reaction involving an ester.
    """
    has_ester_coupling = False

    def dfs_traverse(node):
        nonlocal has_ester_coupling

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ester pattern in reactants
                ester_pattern = Chem.MolFromSmarts("[#6][#8][#6](=[O])")
                amine_pattern = Chem.MolFromSmarts("[N;H2][#6]")

                has_ester = False
                has_amine = False

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if not reactant_mol:
                        continue

                    if reactant_mol.HasSubstructMatch(ester_pattern):
                        has_ester = True

                    if reactant_mol.HasSubstructMatch(amine_pattern):
                        has_amine = True

                # If both ester and amine are present in reactants
                if has_ester and has_amine:
                    # Check if product has fewer fragments than reactants
                    if len(reactants) > 1 and "." not in product:
                        has_ester_coupling = True
                        print(
                            f"Detected fragment coupling with ester at reaction with RSMI: {rsmi}"
                        )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Fragment coupling with ester strategy detected: {has_ester_coupling}")
    return has_ester_coupling

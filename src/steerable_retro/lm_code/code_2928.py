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
    This function detects fragment coupling via isothiocyanate-amine reaction
    to form thiourea.
    """
    has_isothiocyanate_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_isothiocyanate_coupling

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                products = rsmi.split(">")[-1]

                # Check for isothiocyanate-amine coupling
                if "." in reactants:  # Multiple reactants (potential coupling)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants.split(".")]
                    product_mol = Chem.MolFromSmiles(products)

                    if all(reactant_mols) and product_mol:
                        # Isothiocyanate pattern
                        isothiocyanate_pattern = Chem.MolFromSmarts("[#6]-[#7]=[#6]=[#16]")
                        # Amine pattern
                        amine_pattern = Chem.MolFromSmarts("[#7;H2]")
                        # Thiourea pattern
                        thiourea_pattern = Chem.MolFromSmarts("[#7]-[#6](=[#16])-[#7]")

                        has_isothiocyanate = any(
                            mol.HasSubstructMatch(isothiocyanate_pattern) for mol in reactant_mols
                        )
                        has_amine = any(
                            mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols
                        )
                        forms_thiourea = product_mol.HasSubstructMatch(thiourea_pattern)

                        if has_isothiocyanate and has_amine and forms_thiourea:
                            print(f"Detected isothiocyanate-amine coupling at depth {depth}")
                            has_isothiocyanate_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_isothiocyanate_coupling

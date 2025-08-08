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
    Detects if the synthesis involves introduction of a tert-butyl group to an aromatic ring.
    """
    found_tert_butylation = False

    def dfs_traverse(node):
        nonlocal found_tert_butylation

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                product = Chem.MolFromSmiles(product_smiles)

                if product and all(r for r in reactants):
                    # Check for tert-butyl group in product but not in reactants
                    tert_butyl_pattern = Chem.MolFromSmarts("[c]-[C](-[C])(-[C])-[C]")

                    has_tert_butyl_in_product = product.HasSubstructMatch(tert_butyl_pattern)
                    has_tert_butyl_in_reactants = any(
                        r.HasSubstructMatch(tert_butyl_pattern) for r in reactants if r
                    )

                    if has_tert_butyl_in_product and not has_tert_butyl_in_reactants:
                        print(f"Found tert-butylation of aromatic ring")
                        found_tert_butylation = True
            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return found_tert_butylation

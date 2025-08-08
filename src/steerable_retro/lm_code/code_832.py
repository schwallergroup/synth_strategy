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
    This function detects a late-stage ether formation strategy where two complex fragments
    are connected via C-O bond formation in the final steps of the synthesis, with one fragment
    containing a chloromethyl group and the other containing a phenol.
    """
    # Track if we found the strategy
    found_strategy = False

    def dfs_traverse(node, depth=0):
        nonlocal found_strategy

        # Only look at reaction nodes
        if node["type"] == "reaction" and depth <= 1:  # Focus on late-stage reactions (depth 0-1)
            try:
                # Extract reaction SMILES
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if we have at least 2 reactants (convergent step)
                if len(reactants_smiles) >= 2:
                    # Create RDKit mol objects
                    reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Skip if any molecule failed to parse
                    if None in reactant_mols or product_mol is None:
                        print("Failed to parse molecules in reaction")
                        return

                    # Check for chloromethyl group in reactants
                    chloromethyl_pattern = Chem.MolFromSmarts("[Cl][CH2][#6]")
                    has_chloromethyl = any(
                        mol.HasSubstructMatch(chloromethyl_pattern) for mol in reactant_mols
                    )

                    # Check for phenol/alcohol group in reactants
                    phenol_pattern = Chem.MolFromSmarts("[OH][c,C]")
                    has_phenol = any(mol.HasSubstructMatch(phenol_pattern) for mol in reactant_mols)

                    # Check for ether formation in product
                    ether_pattern = Chem.MolFromSmarts("[#6][CH2][O][c,C]")
                    has_ether = product_mol.HasSubstructMatch(ether_pattern)

                    # If we have both required reactant groups and the ether in product, it's our strategy
                    if has_chloromethyl and has_phenol and has_ether:
                        print(f"Found late-stage ether formation at depth {depth}")
                        found_strategy = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_strategy

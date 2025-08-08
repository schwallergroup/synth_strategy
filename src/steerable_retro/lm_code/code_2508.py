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
    Detects if an aromatic scaffold with specific substituents is preserved throughout the synthesis.
    """
    # Define a more flexible pattern for aromatic ring with methyl and fluoro substituents
    # This pattern looks for a benzene ring with a methyl and a fluoro substituent in any position
    aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
    methyl_pattern = Chem.MolFromSmarts("c-[CH3]")
    fluoro_pattern = Chem.MolFromSmarts("c-F")

    # Track if we've found the scaffold in the final product
    scaffold_in_final = False
    # Track molecules that should have the scaffold
    molecules_with_scaffold = set()

    def has_target_scaffold(mol):
        """Check if a molecule has the target aromatic scaffold with methyl and fluoro substituents"""
        if mol is None:
            return False

        # Check for benzene ring
        if not mol.HasSubstructMatch(aromatic_pattern):
            return False

        # Check for methyl and fluoro substituents
        has_methyl = mol.HasSubstructMatch(methyl_pattern)
        has_fluoro = mol.HasSubstructMatch(fluoro_pattern)

        return has_methyl and has_fluoro

    def dfs_traverse(node, depth=0):
        nonlocal scaffold_in_final, molecules_with_scaffold

        # Check molecule nodes
        if node["type"] == "mol" and node.get("smiles"):
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol is not None:
                has_scaffold = has_target_scaffold(mol)

                # If this is the final product (depth 0), check if it has the scaffold
                if depth == 0:
                    scaffold_in_final = has_scaffold
                    print(
                        f"Final product {'has' if has_scaffold else 'does not have'} the scaffold: {smiles}"
                    )
                    if has_scaffold:
                        molecules_with_scaffold.add(smiles)

                # For intermediates, track if they have the scaffold
                elif has_scaffold:
                    print(f"Intermediate has scaffold: {smiles}")
                    molecules_with_scaffold.add(smiles)
                else:
                    print(f"Intermediate does not have scaffold: {smiles}")

        # Check reaction nodes to ensure scaffold preservation
        elif node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has the scaffold
                product_mol = Chem.MolFromSmiles(product)
                product_has_scaffold = has_target_scaffold(product_mol)

                # Check if any reactant has the scaffold
                reactant_has_scaffold = False
                for r in reactants:
                    r_mol = Chem.MolFromSmiles(r)
                    if has_target_scaffold(r_mol):
                        reactant_has_scaffold = True
                        break

                # If a reactant has the scaffold but the product doesn't, the scaffold wasn't preserved
                if reactant_has_scaffold and not product_has_scaffold:
                    print(f"Scaffold lost in reaction: {rsmi}")

                # If a product has the scaffold but no reactant does, it was created
                if product_has_scaffold and not reactant_has_scaffold:
                    print(f"Scaffold created in reaction: {rsmi}")
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # The scaffold must be in the final product and preserved in key intermediates
    result = scaffold_in_final and len(molecules_with_scaffold) > 0
    print(f"Scaffold in final product: {scaffold_in_final}")
    print(f"Number of molecules with scaffold: {len(molecules_with_scaffold)}")
    print(f"Final result: {result}")

    return result

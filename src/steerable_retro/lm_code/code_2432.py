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
    Detects if the synthesis route uses a late-stage amide coupling strategy.
    This checks if an amide bond is formed in the final or penultimate step.
    """
    amide_formation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for amide formation
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(m for m in reactant_mols):
                # Look for amide pattern in product but not in reactants
                amide_pattern = Chem.MolFromSmarts("[C;$(C=O)][N;!$(N=*)]")
                product_matches = product_mol.GetSubstructMatches(amide_pattern)

                # Check if any reactant has the same amide bond
                new_amide_formed = False
                for match in product_matches:
                    amide_bond = (match[0], match[1])
                    # Check if this bond exists in any reactant
                    exists_in_reactants = False
                    for r_mol in reactant_mols:
                        if r_mol.HasSubstructMatch(amide_pattern):
                            r_matches = r_mol.GetSubstructMatches(amide_pattern)
                            for r_match in r_matches:
                                # If same atoms are connected in reactant, it's not a new bond
                                if (
                                    product_mol.GetAtomWithIdx(match[0]).GetSymbol()
                                    == r_mol.GetAtomWithIdx(r_match[0]).GetSymbol()
                                    and product_mol.GetAtomWithIdx(match[1]).GetSymbol()
                                    == r_mol.GetAtomWithIdx(r_match[1]).GetSymbol()
                                ):
                                    exists_in_reactants = True
                                    break

                    if not exists_in_reactants:
                        new_amide_formed = True
                        break

                if new_amide_formed:
                    amide_formation_depth = depth
                    print(f"Amide formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if amide formation occurred in the late stage (depth 0 or 1)
    is_late_stage = amide_formation_depth is not None and amide_formation_depth <= 1
    print(f"Late-stage amide coupling: {is_late_stage} (depth: {amide_formation_depth})")
    return is_late_stage

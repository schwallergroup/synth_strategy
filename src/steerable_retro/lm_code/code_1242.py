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
    This function detects a linear synthesis strategy (as opposed to convergent).

    A linear synthesis is characterized by sequential addition of building blocks,
    where most reaction steps have only one significant reactant (excluding reagents).
    """
    # In a linear synthesis, most reactions have only one non-reagent reactant
    reaction_count = 0
    linear_reaction_count = 0

    # Common reagents to exclude (simplified SMILES)
    common_reagents = {
        "O",
        "CCO",
        "CO",
        "CC(=O)O",
        "C1CCOC1",
        "[Pd]",
        "ClCCl",
        "O=C(O)C(F)(F)F",
        "[Na+]",
        "CC(C)O",
        "CN",
        "CC#N",
        "C1CCCCC1",
        "CS(=O)(=O)O",
        "CC(C)(C)O",
        "O=S(=O)(O)O",
        "N",
        "C(=O)=O",
    }

    def dfs_traverse(node):
        nonlocal reaction_count, linear_reaction_count

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                reaction_count += 1

                # Count significant reactants (excluding small reagents)
                significant_reactants = []
                for r in reactants:
                    # Skip common reagents
                    if r in common_reagents:
                        continue

                    # Skip small molecules (likely reagents)
                    mol = Chem.MolFromSmiles(r)
                    if mol and mol.GetNumHeavyAtoms() > 5:
                        significant_reactants.append(r)

                # Check if this is a linear reaction step
                if len(significant_reactants) <= 1:
                    linear_reaction_count += 1
                    print(f"Linear reaction step detected: {rsmi}")
                else:
                    # Even with multiple significant reactants, check if one is much larger
                    # This could still be a linear step with a significant reagent
                    atom_counts = []
                    for r in significant_reactants:
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            atom_counts.append(mol.GetNumHeavyAtoms())

                    if atom_counts and len(atom_counts) > 1:
                        # If one reactant is at least twice as large as any other,
                        # it's likely the main substrate in a linear synthesis
                        max_atoms = max(atom_counts)
                        others_max = max([c for c in atom_counts if c != max_atoms] or [0])
                        if max_atoms > 2 * others_max:
                            linear_reaction_count += 1
                            print(f"Linear reaction step detected (size disparity): {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If most reactions (>70%) are linear, consider it a linear synthesis
    # For small routes (≤3 reactions), use a lower threshold of 50%
    if reaction_count > 0:
        linear_ratio = linear_reaction_count / reaction_count
        print(f"Linear reaction ratio: {linear_ratio} ({linear_reaction_count}/{reaction_count})")

        if reaction_count <= 3:
            return linear_ratio >= 0.5
        else:
            return linear_ratio >= 0.7

    return False

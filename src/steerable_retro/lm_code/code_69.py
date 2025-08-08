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
    Detects if the synthesis uses a convergent approach with 3+ fragments combined in the final steps.
    A convergent synthesis builds separate complex fragments and combines them in late-stage reactions.
    """
    # Track significant fragments across late-stage reactions
    significant_fragments = set()
    max_depth = 0

    # First, calculate the maximum depth of the synthesis route
    def calculate_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            calculate_max_depth(child, current_depth + 1)

    calculate_max_depth(route)
    print(f"Maximum depth of synthesis route: {max_depth}")

    # Define what constitutes a late-stage reaction (within top 30% of steps)
    late_stage_threshold = max(1, int(max_depth * 0.3))
    print(f"Late-stage threshold depth: {late_stage_threshold}")

    # Helper function to determine if a molecule is a significant fragment
    def is_significant_fragment(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        # Criteria for significant fragments:
        # 1. Must have at least 5 atoms
        # 2. Must have at least 1 carbon atom
        # 3. Must have at least 3 bonds (to filter out simple molecules)
        atom_count = mol.GetNumAtoms()
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        bond_count = mol.GetNumBonds()

        return atom_count >= 5 and carbon_count >= 1 and bond_count >= 3

    # Analyze the synthesis route to find significant fragments in late-stage reactions
    def analyze_route(node, depth=0):
        # Only consider late-stage reactions
        if node["type"] == "reaction" and depth <= late_stage_threshold:
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                # Add significant fragments to our set
                for reactant in reactants:
                    if is_significant_fragment(reactant):
                        # Use canonical SMILES to avoid duplicates due to different representations
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                            significant_fragments.add(canonical_smiles)

                print(
                    f"Depth {depth}: Found {len(reactants)} reactants in reaction, "
                    f"{sum(1 for r in reactants if is_significant_fragment(r))} significant"
                )

        # Continue traversal
        for child in node.get("children", []):
            analyze_route(child, depth + 1)

    # Start analysis
    analyze_route(route)

    # Check if we have 3 or more unique significant fragments
    print(
        f"Total unique significant fragments found in late-stage reactions: {len(significant_fragments)}"
    )

    is_convergent = len(significant_fragments) >= 3
    if is_convergent:
        print("Found convergent synthesis strategy with 3+ significant fragments")
    else:
        print(
            "Not a convergent synthesis: fewer than 3 significant fragments in late-stage reactions"
        )

    return is_convergent

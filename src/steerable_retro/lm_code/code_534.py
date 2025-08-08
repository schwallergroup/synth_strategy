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
    This function detects if the synthesis follows a linear strategy without convergent steps.
    A linear synthesis has:
    1. Each reaction step has at most 2 significant reactants
    2. The product of each step becomes a reactant in the next step
    3. At least 3 reaction steps in total
    """
    is_linear = True
    reaction_nodes = []

    # First pass: collect all reaction nodes in order
    def collect_reactions(node, depth=0):
        if node["type"] == "reaction":
            reaction_nodes.append((node, depth))

        for child in node.get("children", []):
            collect_reactions(child, depth + 1)

    collect_reactions(route)

    # Sort reactions by depth (from late-stage to early-stage)
    reaction_nodes.sort(key=lambda x: x[1])

    print(f"Total reaction nodes found: {len(reaction_nodes)}")

    # Check if any reaction has more than 2 significant reactants
    for node, depth in reaction_nodes:
        if "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]
            reactants = reactants_part.split(".")

            # Filter out small molecules (likely reagents/solvents)
            significant_reactants = []
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.GetNumHeavyAtoms() > 3:  # Consider only significant molecules
                    significant_reactants.append(reactant)

            print(
                f"Reaction at depth {depth} has {len(significant_reactants)} significant reactants"
            )

            # If any step has more than 2 significant reactants, it's not a linear synthesis
            if len(significant_reactants) > 2:
                print(
                    f"Non-linear step detected with {len(significant_reactants)} significant reactants"
                )
                is_linear = False
                break

    # Check if products flow into subsequent reactions (true linear path)
    if is_linear and len(reaction_nodes) >= 2:
        for i in range(len(reaction_nodes) - 1):
            current_node, current_depth = reaction_nodes[i]
            next_node, next_depth = reaction_nodes[i + 1]

            # Get current product
            current_rsmi = current_node["metadata"]["rsmi"]
            current_product = current_rsmi.split(">")[-1]
            current_mol = Chem.MolFromSmiles(current_product)

            # Get next reactants
            next_rsmi = next_node["metadata"]["rsmi"]
            next_reactants = next_rsmi.split(">")[0].split(".")

            # Check if current product is used in next reaction
            product_used = False

            for reactant in next_reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and current_mol:
                    # Use MCS to compare structures with a threshold
                    if reactant_mol.GetNumHeavyAtoms() > 3:  # Only consider significant reactants
                        # Try exact match first (faster)
                        if Chem.MolToSmiles(current_mol, isomericSmiles=False) == Chem.MolToSmiles(
                            reactant_mol, isomericSmiles=False
                        ):
                            product_used = True
                            print(
                                f"Product from depth {current_depth} matches reactant in depth {next_depth} (exact match)"
                            )
                            break

                        # If exact match fails, try MCS
                        mcs = rdFMCS.FindMCS(
                            [current_mol, reactant_mol],
                            atomCompare=rdFMCS.AtomCompare.CompareElements,
                            bondCompare=rdFMCS.BondCompare.CompareOrder,
                            completeRingsOnly=True,
                            ringMatchesRingOnly=True,
                            matchValences=False,
                        )

                        if (
                            mcs.numAtoms
                            >= min(current_mol.GetNumHeavyAtoms(), reactant_mol.GetNumHeavyAtoms())
                            * 0.8
                        ):
                            product_used = True
                            print(
                                f"Product from depth {current_depth} matches reactant in depth {next_depth} (MCS match: {mcs.numAtoms} atoms)"
                            )
                            break

            if not product_used:
                print(
                    f"Break in linear path detected: product from depth {current_depth} not used in reaction at depth {next_depth}"
                )
                is_linear = False
                break

    print(f"Reaction count: {len(reaction_nodes)}")
    print(f"Linear synthesis strategy: {is_linear}")

    # A true linear synthesis should have multiple steps
    return is_linear and len(reaction_nodes) >= 3

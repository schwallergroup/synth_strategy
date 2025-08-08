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
    This function detects a linear synthesis strategy where each step has only one product
    that serves as a reactant in the next step, as opposed to convergent synthesis.

    In a linear synthesis:
    1. Each intermediate molecule is used in exactly one reaction
    2. The synthesis follows a sequential path without branching/converging
    """
    # Track molecules and their reaction parents
    molecule_parents = {}
    # Track reaction nodes for counting
    reaction_nodes = []
    # Track the target molecule (root of the tree)
    target_molecule = None

    def dfs_traverse(node, parent=None):
        nonlocal target_molecule

        # If this is the root node and it's a molecule, it's our target
        if parent is None and node["type"] == "mol":
            target_molecule = node["smiles"]

        if node["type"] == "mol":
            # For molecule nodes, track their reaction parent
            mol_smiles = node["smiles"]
            if parent is not None and parent["type"] == "reaction":
                if mol_smiles in molecule_parents:
                    molecule_parents[mol_smiles].append(parent)
                else:
                    molecule_parents[mol_smiles] = [parent]

        elif node["type"] == "reaction":
            # Track reaction nodes
            reaction_nodes.append(node)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, node)

    # Start traversal
    dfs_traverse(route)

    # If there's only one or zero reactions, it's linear by definition
    if len(reaction_nodes) <= 1:
        return True

    # Check if each molecule (except the target and starting materials) has exactly one reaction parent
    is_linear = True
    for mol_smiles, parents in molecule_parents.items():
        # Skip the target molecule and starting materials
        if mol_smiles == target_molecule or any(
            node.get("in_stock", False)
            for node in route.get("children", [])
            if node["type"] == "mol" and node["smiles"] == mol_smiles
        ):
            continue

        # If a molecule has more than one reaction parent, the synthesis is not linear
        if len(parents) > 1:
            print(
                f"Molecule {mol_smiles} is used in multiple reactions, indicating convergent synthesis"
            )
            is_linear = False
            break

    # Additional check: ensure reactions form a single path
    # Sort reactions by their position in the synthesis
    if is_linear and len(reaction_nodes) > 1:
        # Build a directed graph of reactions
        reaction_graph = {}
        for i, rxn in enumerate(reaction_nodes):
            reaction_graph[i] = []

            # Get the product of this reaction
            try:
                rsmi = rxn["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Find reactions that use this product
                for j, other_rxn in enumerate(reaction_nodes):
                    if i == j:
                        continue

                    try:
                        other_rsmi = other_rxn["metadata"]["rsmi"]
                        reactants = other_rsmi.split(">")[0].split(".")

                        if product in reactants:
                            reaction_graph[i].append(j)
                    except:
                        pass
            except:
                pass

        # Check if the graph forms a single path
        # In a linear synthesis, each reaction should have at most one outgoing edge
        for rxn_idx, next_rxns in reaction_graph.items():
            if len(next_rxns) > 1:
                print(
                    f"Reaction {rxn_idx} leads to multiple subsequent reactions, indicating branching"
                )
                is_linear = False
                break

    return is_linear

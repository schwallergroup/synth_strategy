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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    This function detects a convergent synthesis strategy where two heterocycle-forming
    pathways converge before a final bond-forming step.
    """
    # Initialize tracking variables
    heterocycle_formations = []  # (depth, branch, ring, mol_smiles)
    convergent_reactions = []  # (depth, reaction_node, [branch_ids])
    branch_paths = {}  # branch_id -> list of molecule smiles in the path

    # List of heterocyclic rings to check
    heterocycles = [
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "quinoline",
        "isoquinoline",
        "indole",
        "purine",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "furan",
        "thiophene",
        "isoxazole",
        "isothiazole",
    ]

    # Track molecule nodes to their branch IDs
    mol_to_branch = {}
    next_branch_id = 0

    def get_branch_id(mol_smiles):
        nonlocal next_branch_id, mol_to_branch
        if mol_smiles not in mol_to_branch:
            mol_to_branch[mol_smiles] = next_branch_id
            branch_paths[next_branch_id] = [mol_smiles]
            next_branch_id += 1
        return mol_to_branch[mol_smiles]

    def dfs_traverse(node, depth=0, parent_branch=None):
        nonlocal heterocycle_formations, convergent_reactions

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            current_branch = get_branch_id(mol_smiles)

            # If this is a continuation of a branch, update the branch path
            if parent_branch is not None and parent_branch == current_branch:
                if mol_smiles not in branch_paths[current_branch]:
                    branch_paths[current_branch].append(mol_smiles)

            # Check for heterocycles in this molecule
            mol_heterocycles = []
            for ring in heterocycles:
                if checker.check_ring(ring, mol_smiles):
                    mol_heterocycles.append(ring)

            # Traverse children with the current branch ID
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1, current_branch)

        elif node["type"] == "reaction":
            try:
                # Extract reaction information
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Get branch IDs of reactants
                reactant_branches = []
                for r_smiles in reactants_smiles:
                    if r_smiles in mol_to_branch:
                        reactant_branches.append(mol_to_branch[r_smiles])
                    else:
                        # If reactant not seen before, assign a new branch ID
                        branch_id = get_branch_id(r_smiles)
                        reactant_branches.append(branch_id)

                # Check for heterocycle formation
                for ring in heterocycles:
                    # Check if product contains heterocycle
                    if checker.check_ring(ring, product_smiles):
                        # Check which reactants don't have the heterocycle
                        for i, r_smiles in enumerate(reactants_smiles):
                            if not checker.check_ring(ring, r_smiles):
                                branch_id = (
                                    reactant_branches[i]
                                    if i < len(reactant_branches)
                                    else get_branch_id(r_smiles)
                                )
                                print(
                                    f"Detected {ring} formation at depth {depth} in branch {branch_id}"
                                )
                                heterocycle_formations.append(
                                    (depth, branch_id, ring, product_smiles)
                                )

                # Check for convergent pattern - reaction with reactants from different branches
                unique_branches = set(reactant_branches)
                if len(unique_branches) >= 2:
                    print(
                        f"Detected convergent reaction at depth {depth} with branches {unique_branches}"
                    )
                    convergent_reactions.append((depth, node, list(unique_branches)))
            except Exception as e:
                print(f"Error processing reaction node: {e}")

            # Traverse children
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1, parent_branch)

    # Start traversal from root
    dfs_traverse(route)

    # Get branches that form heterocycles
    heterocycle_branches = set(branch for _, branch, _, _ in heterocycle_formations)
    print(f"Branches with heterocycle formation: {heterocycle_branches}")

    # Check if we have a convergent reaction that combines heterocycle-forming branches
    convergent_heterocycle_synthesis = False
    for depth, _, branches in convergent_reactions:
        # Check if this convergent reaction combines at least two heterocycle-forming branches
        heterocycle_branches_in_reaction = [b for b in branches if b in heterocycle_branches]
        if len(heterocycle_branches_in_reaction) >= 2:
            print(
                f"Found convergent reaction at depth {depth} combining heterocycle branches: {heterocycle_branches_in_reaction}"
            )
            convergent_heterocycle_synthesis = True
            break

    # If no direct convergent reaction found, check if heterocycle-forming branches eventually converge
    if not convergent_heterocycle_synthesis and len(heterocycle_branches) >= 2:
        # Check if any convergent reaction has ancestors from different heterocycle-forming branches
        for depth, node, branches in convergent_reactions:
            # Check if branches in this convergent reaction have ancestors from heterocycle-forming branches
            branch_ancestors = set()
            for branch in branches:
                # Check if this branch or any of its ancestors formed a heterocycle
                if branch in heterocycle_branches:
                    branch_ancestors.add(branch)
                # Check branch path for molecules that are in heterocycle-forming branches
                for mol_smiles in branch_paths.get(branch, []):
                    for hc_depth, hc_branch, _, hc_mol in heterocycle_formations:
                        if mol_smiles == hc_mol and hc_branch != branch:
                            branch_ancestors.add(hc_branch)

            if len(branch_ancestors.intersection(heterocycle_branches)) >= 2:
                print(
                    f"Found convergent reaction at depth {depth} with ancestors from heterocycle branches: {branch_ancestors}"
                )
                convergent_heterocycle_synthesis = True
                break

    print(f"Heterocycle formations: {heterocycle_formations}")
    print(f"Convergent reactions: {convergent_reactions}")
    print(f"Strategy detected: {convergent_heterocycle_synthesis}")

    return convergent_heterocycle_synthesis

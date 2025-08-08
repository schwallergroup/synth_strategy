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
    This function detects a strategy where a nitro-substituted aromatic ring
    remains unchanged throughout the synthesis while other functional groups
    are modified.
    """
    # Track main synthetic pathway molecules
    main_pathway_mols = []
    reaction_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count

        if node["type"] == "mol":
            # Only track the main product (target molecule or intermediate)
            if depth == 0 or node.get("in_stock", False) == False:
                main_pathway_mols.append(node["smiles"])
                print(f"Added to main pathway: {node['smiles']}")

        elif node["type"] == "reaction":
            reaction_count += 1

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if all main pathway molecules have nitro on aromatic ring
    if len(main_pathway_mols) < 2:
        print("Not enough molecules in main pathway")
        return False

    # Check if all molecules have nitro group on aromatic ring
    all_have_nitro_aromatic = True
    for mol_smiles in main_pathway_mols:
        has_nitro = checker.check_fg("Nitro group", mol_smiles)

        # Check if the nitro group is on an aromatic ring
        mol = Chem.MolFromSmiles(mol_smiles)
        nitro_on_aromatic = False

        if has_nitro and mol:
            print(f"Checking nitro on aromatic for: {mol_smiles}")

            # Get all aromatic atoms
            aromatic_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()]

            # For each aromatic atom, check if it's connected to a nitro group
            for aromatic_idx in aromatic_atoms:
                aromatic_atom = mol.GetAtomWithIdx(aromatic_idx)
                for neighbor in aromatic_atom.GetNeighbors():
                    # If the neighbor is a nitrogen with a double bond to oxygen, it's likely part of a nitro group
                    if neighbor.GetSymbol() == "N":
                        for n_neighbor in neighbor.GetNeighbors():
                            if (
                                n_neighbor.GetSymbol() == "O"
                                and n_neighbor.GetIdx() != aromatic_idx
                            ):
                                nitro_on_aromatic = True
                                break
                        if nitro_on_aromatic:
                            break
                if nitro_on_aromatic:
                    break

        if not (has_nitro and nitro_on_aromatic):
            all_have_nitro_aromatic = False
            print(f"Molecule without nitro on aromatic ring: {mol_smiles}")
            break

    # Strategy is detected if all main pathway molecules have nitro on aromatic ring
    # and there are at least 3 reaction steps
    strategy_detected = all_have_nitro_aromatic and reaction_count >= 3
    print(f"Persistent nitro aromatic strategy detected: {strategy_detected}")
    print(f"Reaction count: {reaction_count}")
    print(f"Main pathway molecules: {len(main_pathway_mols)}")

    return strategy_detected

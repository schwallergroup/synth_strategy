#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects if C-C bond formation occurs in the late stage of synthesis (low depth).
    """
    cc_bond_formation_depths = []
    max_depth = 0

    # List of C-C bond forming reactions
    cc_bond_reactions = [
        "Suzuki",
        "Negishi",
        "Heck",
        "Stille",
        "Sonogashira",
        "Kumada",
        "Hiyama-Denmark",
        "Friedel-Crafts alkylation",
        "Diels-Alder",
        "Wittig",
        "Grignard",
        "Aldol condensation",
        "Michael addition",
        "Aryllithium cross-coupling",
        "Catellani",
        "decarboxylative_coupling",
        "A3 coupling",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a C-C bond forming reaction
            is_cc_formation = False

            # Method 1: Check against known C-C bond forming reactions
            for rxn_type in cc_bond_reactions:
                try:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Found {rxn_type} reaction at depth {depth}: {rsmi}")
                        is_cc_formation = True
                        break
                except Exception as e:
                    print(f"Error checking reaction type {rxn_type}: {e}")

            # Method 2: If not identified by reaction type, check for new C-C bonds using atom mapping
            if not is_cc_formation:
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                    if product_mol and all(reactant_mols):
                        # Extract atom maps from reactants
                        reactant_c_maps = {}
                        for r_mol in reactant_mols:
                            for atom in r_mol.GetAtoms():
                                if (
                                    atom.GetAtomicNum() == 6
                                    and atom.GetAtomMapNum() > 0
                                ):
                                    reactant_c_maps[
                                        atom.GetAtomMapNum()
                                    ] = atom.GetIdx()

                        # Check for new C-C bonds in product between mapped atoms
                        for bond in product_mol.GetBonds():
                            begin_atom = bond.GetBeginAtom()
                            end_atom = bond.GetEndAtom()

                            if (
                                begin_atom.GetAtomicNum() == 6
                                and end_atom.GetAtomicNum() == 6
                                and begin_atom.GetAtomMapNum() > 0
                                and end_atom.GetAtomMapNum() > 0
                            ):

                                map1 = begin_atom.GetAtomMapNum()
                                map2 = end_atom.GetAtomMapNum()

                                # Check if these mapped atoms were in different reactant molecules
                                # or if this bond didn't exist in the reactants
                                new_bond = True
                                for r_mol in reactant_mols:
                                    atom1_idx = None
                                    atom2_idx = None

                                    for atom in r_mol.GetAtoms():
                                        if atom.GetAtomMapNum() == map1:
                                            atom1_idx = atom.GetIdx()
                                        elif atom.GetAtomMapNum() == map2:
                                            atom2_idx = atom.GetIdx()

                                    if atom1_idx is not None and atom2_idx is not None:
                                        bond = r_mol.GetBondBetweenAtoms(
                                            atom1_idx, atom2_idx
                                        )
                                        if bond is not None:
                                            new_bond = False
                                            break

                                if new_bond:
                                    print(
                                        f"Found new C-C bond formation at depth {depth} between atoms with maps {map1} and {map2}"
                                    )
                                    is_cc_formation = True
                                    break
                except Exception as e:
                    print(f"Error analyzing atom mapping: {e}")

            if is_cc_formation:
                cc_bond_formation_depths.append(depth)

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Define late stage as the first third of the synthesis depth
    late_stage_threshold = max(1, max_depth // 3)
    result = any(depth <= late_stage_threshold for depth in cc_bond_formation_depths)

    print(f"Max depth: {max_depth}, Late stage threshold: {late_stage_threshold}")
    print(f"C-C bond formation depths: {cc_bond_formation_depths}")
    print(f"Late-stage C-C bond formation: {result}")

    return result

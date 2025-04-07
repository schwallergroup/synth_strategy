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
    This function detects the preservation of a dichlorobenzyl moiety throughout the synthesis.
    """

    def has_dichlorobenzyl_moiety(smiles):
        """Check if a molecule contains a dichlorobenzyl moiety"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        # Check for aromatic ring with two chlorines and a CH2 group
        # First check if molecule has aromatic halides
        if checker.check_fg("Aromatic halide", smiles):
            # Get the indices of aromatic rings
            ring_info = mol.GetRingInfo()
            if ring_info.NumRings() > 0:
                # Check each aromatic ring for two chlorines
                for ring_atoms in ring_info.AtomRings():
                    # Check if ring is aromatic
                    if all(
                        mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms
                    ):
                        # Count chlorines attached to the ring
                        chlorine_count = 0
                        benzyl_found = False

                        for atom_idx in ring_atoms:
                            atom = mol.GetAtomWithIdx(atom_idx)

                            # Check neighbors for chlorines
                            for neighbor in atom.GetNeighbors():
                                if neighbor.GetSymbol() == "Cl":
                                    chlorine_count += 1

                                # Check for CH2 group (benzyl)
                                if (
                                    neighbor.GetSymbol() == "C"
                                    and not neighbor.GetIsAromatic()
                                ):
                                    # Verify it's a CH2 group
                                    if neighbor.GetTotalNumHs() == 2 or (
                                        neighbor.GetTotalNumHs()
                                        + neighbor.GetNumExplicitHs()
                                        == 2
                                    ):
                                        benzyl_found = True

                        if chlorine_count >= 2 and benzyl_found:
                            return True

        return False

    # Check if the target molecule contains a dichlorobenzyl moiety
    if route["type"] == "mol" and not has_dichlorobenzyl_moiety(route["smiles"]):
        print(f"Target molecule does not contain dichlorobenzyl moiety")
        return False

    # Track if all reactions preserve the moiety
    all_preserved = True
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal all_preserved, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product has dichlorobenzyl moiety
                product_has_moiety = has_dichlorobenzyl_moiety(product_smiles)

                # Check if at least one reactant has dichlorobenzyl moiety
                reactant_has_moiety = any(
                    has_dichlorobenzyl_moiety(r) for r in reactants_smiles
                )

                # If moiety is in product but not in any reactant, or vice versa, it's not preserved
                if product_has_moiety != reactant_has_moiety:
                    print(f"Dichlorobenzyl moiety not preserved in reaction: {rsmi}")
                    all_preserved = False
                elif product_has_moiety:
                    print(f"Dichlorobenzyl moiety preserved in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return true if all reactions preserve the moiety and at least one reaction exists
    return all_preserved and reaction_count > 0

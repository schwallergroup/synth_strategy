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
    This function detects if the synthesis includes C-C bond formation via
    aryl-alkyne coupling using a triflate leaving group.
    """
    coupling_found = False

    def dfs_traverse(node):
        nonlocal coupling_found

        if node["type"] == "reaction" and not coupling_found:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for triflate-based Sonogashira coupling
                if checker.check_reaction("Sonogashira alkyne_aryl OTf", rsmi):
                    print(f"Found aryl-alkyne coupling using triflate: {rsmi}")
                    coupling_found = True
                    return

                # If specific reaction check fails, check for components
                triflate_present = any(
                    checker.check_fg("Triflate", reactant) for reactant in reactants
                )
                alkyne_present = any(
                    checker.check_fg("Alkyne", reactant) for reactant in reactants
                )

                # Check if product contains new C-C bond between aryl and alkyne
                if triflate_present and alkyne_present:
                    # Check if any reactant has triflate on aromatic ring
                    aryl_triflate = False
                    for reactant in reactants:
                        if checker.check_fg("Triflate", reactant):
                            # Check if triflate is attached to aromatic ring
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                for atom in mol.GetAtoms():
                                    if atom.GetIsAromatic() and any(
                                        nbr.GetSymbol() == "O"
                                        and any(
                                            nnbr.GetSymbol() == "S"
                                            for nnbr in nbr.GetNeighbors()
                                        )
                                        for nbr in atom.GetNeighbors()
                                    ):
                                        aryl_triflate = True
                                        break

                    # Check if product has aryl-alkyne bond
                    if aryl_triflate:
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            for bond in product_mol.GetBonds():
                                a1 = bond.GetBeginAtom()
                                a2 = bond.GetEndAtom()
                                # Check if bond is between aromatic atom and alkyne carbon
                                if (
                                    a1.GetIsAromatic()
                                    and any(
                                        b.GetBondType() == Chem.BondType.TRIPLE
                                        for b in a2.GetBonds()
                                    )
                                ) or (
                                    a2.GetIsAromatic()
                                    and any(
                                        b.GetBondType() == Chem.BondType.TRIPLE
                                        for b in a1.GetBonds()
                                    )
                                ):
                                    print(
                                        f"Found aryl-alkyne coupling using triflate (component analysis): {rsmi}"
                                    )
                                    coupling_found = True
                                    return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return coupling_found

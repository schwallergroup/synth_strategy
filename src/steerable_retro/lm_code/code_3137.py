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
    This function detects if the synthetic route involves biphenyl formation via cross-coupling
    early in the synthesis (higher depth values).
    """
    cross_coupling_detected = False

    def dfs_traverse(node, current_depth=0):
        nonlocal cross_coupling_detected

        # Add depth to metadata for future reference
        if "metadata" not in node:
            node["metadata"] = {}
        node["metadata"]["depth"] = current_depth

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a cross-coupling reaction
            is_suzuki = (
                checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
            )
            is_negishi = checker.check_reaction("Negishi coupling", rsmi)
            is_stille = any(
                checker.check_reaction(f"Stille reaction_{t}", rsmi)
                for t in [
                    "aryl",
                    "aryl OTf",
                    "vinyl",
                    "vinyl OTf",
                    "benzyl",
                    "benzyl OTf",
                    "allyl",
                    "allyl OTf",
                ]
            )
            is_hiyama = checker.check_reaction("Hiyama-Denmark Coupling", rsmi)
            is_kumada = checker.check_reaction("Kumada cross-coupling", rsmi)
            is_ullmann = checker.check_reaction("Ullmann condensation", rsmi)

            is_cross_coupling = (
                is_suzuki or is_negishi or is_stille or is_hiyama or is_kumada or is_ullmann
            )

            # If not detected by reaction checkers, try to detect by functional groups
            if not is_cross_coupling:
                # Look for aryl halide in one reactant
                has_aryl_halide = any(checker.check_fg("Aromatic halide", r) for r in reactants)
                # Look for boronic acid/ester in another reactant
                has_boronic = any(
                    checker.check_fg("Boronic acid", r) or checker.check_fg("Boronic ester", r)
                    for r in reactants
                )
                # Look for organometallic reagents
                has_organometallic = any(
                    checker.check_fg("Magnesium halide", r)
                    or checker.check_fg("Zinc halide", r)
                    or checker.check_fg("Tin", r)
                    for r in reactants
                )

                is_cross_coupling = has_aryl_halide and (has_boronic or has_organometallic)

            # Check if product has biphenyl structure
            if is_cross_coupling:
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    # Specific biphenyl pattern: two benzene rings directly connected
                    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c1ccccc1")

                    # Check if biphenyl exists in product
                    if prod_mol and prod_mol.HasSubstructMatch(biphenyl_pattern):
                        # Check if biphenyl is formed in the reaction (not pre-existing)
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                        biphenyl_in_reactants = any(
                            mol and mol.HasSubstructMatch(biphenyl_pattern)
                            for mol in reactant_mols
                            if mol
                        )

                        if not biphenyl_in_reactants:
                            # Check if this is early in the synthesis (depth >= 2)
                            depth = current_depth
                            if depth >= 2:
                                cross_coupling_detected = True
                                print(
                                    f"Detected biphenyl formation via cross-coupling at depth {depth}"
                                )
                                print(f"Reaction SMILES: {rsmi}")
                except Exception as e:
                    print(f"Error processing molecule: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    dfs_traverse(route)
    return cross_coupling_detected

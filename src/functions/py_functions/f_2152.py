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
    Detects lactam ring formation via C-N bond formation.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Skip if already found a match
            if result:
                return

            # Check if this is a ring formation reaction
            reactants_mol = Chem.MolFromSmiles(reactants_part)
            product_mol = Chem.MolFromSmiles(product_part)

            if reactants_mol and product_mol:
                # Check if product has more rings than reactants
                reactant_rings = reactants_mol.GetRingInfo().NumRings()
                product_rings = product_mol.GetRingInfo().NumRings()

                if product_rings > reactant_rings:
                    # Check for lactam rings in the product
                    lactam_found = False

                    # Check for common lactam-containing rings
                    if checker.check_ring("pyrrolidone", product_part):
                        print("Found pyrrolidone ring in product")
                        lactam_found = True
                    elif checker.check_ring(
                        "piperidine", product_part
                    ) and checker.check_fg("Tertiary amide", product_part):
                        print("Found piperidine ring with amide in product")
                        lactam_found = True

                    # If no specific lactam ring found, check for amide in ring
                    if not lactam_found:
                        # Check for amide functional group in a ring
                        if (
                            checker.check_fg("Primary amide", product_part)
                            or checker.check_fg("Secondary amide", product_part)
                            or checker.check_fg("Tertiary amide", product_part)
                        ):

                            # Verify the amide is in a ring
                            amide_pattern = Chem.MolFromSmarts("C(=O)N")
                            if product_mol.HasSubstructMatch(amide_pattern):
                                matches = product_mol.GetSubstructMatches(amide_pattern)
                                for match in matches:
                                    # Check if these atoms are in a ring
                                    ring_info = product_mol.GetRingInfo()
                                    if all(
                                        ring_info.IsAtomInRingOfSize(atom_idx, 0)
                                        for atom_idx in match
                                    ):
                                        print("Found amide in a ring structure")
                                        lactam_found = True
                                        break

                    # Verify the amide wasn't present in reactants
                    if lactam_found:
                        reactant_has_lactam = False
                        if (
                            checker.check_fg("Primary amide", reactants_part)
                            or checker.check_fg("Secondary amide", reactants_part)
                            or checker.check_fg("Tertiary amide", reactants_part)
                        ):

                            # Check if the amide is in a ring in reactants
                            amide_pattern = Chem.MolFromSmarts("C(=O)N")
                            if reactants_mol.HasSubstructMatch(amide_pattern):
                                matches = reactants_mol.GetSubstructMatches(
                                    amide_pattern
                                )
                                for match in matches:
                                    ring_info = reactants_mol.GetRingInfo()
                                    if all(
                                        ring_info.IsAtomInRingOfSize(atom_idx, 0)
                                        for atom_idx in match
                                    ):
                                        reactant_has_lactam = True
                                        break

                        # Check for reaction types consistent with lactam formation
                        is_lactam_forming_reaction = False
                        if (
                            checker.check_reaction(
                                "Intramolecular transesterification/Lactone formation",
                                rsmi,
                            )
                            or checker.check_reaction(
                                "Formation of NOS Heterocycles", rsmi
                            )
                            or checker.check_reaction(
                                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                                rsmi,
                            )
                        ):
                            is_lactam_forming_reaction = True

                        # If lactam is new and reaction type is consistent, we found a match
                        if not reactant_has_lactam and (
                            is_lactam_forming_reaction
                            or not checker.check_reaction(
                                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                                rsmi,
                            )
                        ):
                            print("Found lactam ring formation via C-N bond formation")
                            result = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return result

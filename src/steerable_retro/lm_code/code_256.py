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
    Detects if the synthesis uses a Suzuki coupling (aryl halide + boronic acid)
    as the final step (depth 0) to form a biaryl system.
    """
    final_step_is_suzuki = False

    def dfs_traverse(node):
        nonlocal final_step_is_suzuki

        if node["type"] == "reaction":
            # Calculate if this is a final step (all children are molecules)
            is_final_step = node.get("children", []) and all(
                child.get("type") == "mol" for child in node.get("children", [])
            )
            # Also check metadata depth if available
            metadata_depth = node.get("metadata", {}).get("depth", -1)

            print(
                f"Examining reaction node. Calculated final step: {is_final_step}, Metadata depth: {metadata_depth}"
            )

            if is_final_step or metadata_depth == 0:
                # Get reaction SMILES
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    print("No reaction SMILES found in metadata")
                    return

                print(f"Reaction SMILES: {rsmi}")

                # Check if this is a Suzuki coupling using the checker function
                is_suzuki_acids = checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                is_suzuki_esters = checker.check_reaction(
                    "Suzuki coupling with boronic esters", rsmi
                )
                is_suzuki_otf_acids = checker.check_reaction(
                    "Suzuki coupling with boronic acids OTf", rsmi
                )
                is_suzuki_otf_esters = checker.check_reaction(
                    "Suzuki coupling with boronic esters OTf", rsmi
                )
                is_suzuki_sulfonic = checker.check_reaction(
                    "Suzuki coupling with sulfonic esters", rsmi
                )

                print(
                    f"Suzuki checks - acids: {is_suzuki_acids}, esters: {is_suzuki_esters}, otf_acids: {is_suzuki_otf_acids}, otf_esters: {is_suzuki_otf_esters}, sulfonic: {is_suzuki_sulfonic}"
                )

                if any(
                    [
                        is_suzuki_acids,
                        is_suzuki_esters,
                        is_suzuki_otf_acids,
                        is_suzuki_otf_esters,
                        is_suzuki_sulfonic,
                    ]
                ):
                    print("Detected late-stage Suzuki coupling")
                    final_step_is_suzuki = True
                else:
                    print("Not a recognized Suzuki coupling reaction, performing additional checks")

                    # Additional check for biaryl formation
                    try:
                        reactants_part = rsmi.split(">")[0]
                        products_part = rsmi.split(">")[-1]

                        # Check for boronic acid/ester in reactants
                        reactants = reactants_part.split(".")
                        has_boronic = False
                        has_aryl_halide = False

                        for reactant in reactants:
                            if checker.check_fg("Boronic acid", reactant) or checker.check_fg(
                                "Boronic ester", reactant
                            ):
                                print(f"Found boronic acid/ester in reactant: {reactant}")
                                has_boronic = True
                            if checker.check_fg("Aromatic halide", reactant) or checker.check_fg(
                                "Triflate", reactant
                            ):
                                print(f"Found aryl halide/triflate in reactant: {reactant}")
                                has_aryl_halide = True

                        # Check for biaryl in product
                        product_mol = Chem.MolFromSmiles(products_part)
                        if product_mol:
                            biaryl_pattern = Chem.MolFromSmarts("c:c-c:c")
                            has_biaryl = product_mol.HasSubstructMatch(biaryl_pattern)
                            print(f"Product has biaryl: {has_biaryl}")

                            if has_boronic and has_aryl_halide and has_biaryl:
                                print(
                                    "Detected potential Suzuki coupling based on reactants and product structure"
                                )
                                final_step_is_suzuki = True
                    except Exception as e:
                        print(f"Error in additional check: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Final result: {final_step_is_suzuki}")
    return final_step_is_suzuki

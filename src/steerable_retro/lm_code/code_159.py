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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects if the synthetic route involves O-alkylation of a piperidine ring.
    """
    o_alkylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal o_alkylation_detected

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if product contains piperidine ring
            if checker.check_ring("piperidine", product):
                print(f"Piperidine ring found in product at depth {depth}")

                # Check if this is an O-alkylation reaction
                is_o_alkylation_reaction = (
                    checker.check_reaction("Williamson Ether Synthesis", rsmi)
                    or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                    or checker.check_reaction("Mitsunobu_phenole", rsmi)
                    or checker.check_reaction("Williamson ether", rsmi)
                    or checker.check_reaction("Mitsunobu_imide", rsmi)
                    or checker.check_reaction("Williamson Ether Synthesis (intra to epoxy)", rsmi)
                )

                if is_o_alkylation_reaction:
                    print(f"O-alkylation reaction detected at depth {depth}")

                    # Check if piperidine was in reactants
                    piperidine_in_reactants = any(
                        checker.check_ring("piperidine", r) for r in reactants
                    )

                    if piperidine_in_reactants:
                        print(f"Piperidine was present in reactants at depth {depth}")
                        o_alkylation_detected = True
                        return

                # Check for ether formation
                if checker.check_fg("Ether", product):
                    print(f"Ether found in product at depth {depth}")

                    # Check if piperidine was in reactants
                    piperidine_in_reactants = False
                    alcohol_in_reactants = False

                    for reactant in reactants:
                        if checker.check_ring("piperidine", reactant):
                            piperidine_in_reactants = True
                            # Check if this piperidine has an alcohol group
                            if (
                                checker.check_fg("Primary alcohol", reactant)
                                or checker.check_fg("Secondary alcohol", reactant)
                                or checker.check_fg("Tertiary alcohol", reactant)
                                or checker.check_fg("Aromatic alcohol", reactant)
                                or checker.check_fg("Phenol", reactant)
                            ):
                                alcohol_in_reactants = True
                                print(
                                    f"Piperidine with alcohol found in reactants at depth {depth}"
                                )
                                break

                    if piperidine_in_reactants and alcohol_in_reactants:
                        # Check if the alcohol was converted to an ether
                        # This is a key step in O-alkylation
                        print(f"Potential piperidine O-alkylation detected at depth {depth}")

                        # Additional check for common O-alkylation reagents
                        reagents = rsmi.split(">")[1].split(".")
                        mitsunobu_reagents = any(
                            "P(" in r or "DEAD" in r or "DIAD" in r for r in reagents
                        )

                        if mitsunobu_reagents or is_o_alkylation_reaction:
                            print(f"Confirmed piperidine O-alkylation at depth {depth}")
                            o_alkylation_detected = True
                            return

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: Piperidine O-alkylation detected = {o_alkylation_detected}")

    return o_alkylation_detected

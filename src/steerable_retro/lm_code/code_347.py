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
    This function detects if the synthetic route involves a late-stage cross-coupling reaction.
    Late-stage is defined as the final reaction step (depth 0) or the step immediately before it (depth 1).
    """
    has_late_stage_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_coupling

        if node["type"] == "reaction" and depth <= 1:  # Late stage = final step or step before
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Check for cross-coupling reactions directly
                if (
                    checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                    or checker.check_reaction("Negishi coupling", rsmi)
                    or checker.check_reaction("Stille reaction_aryl", rsmi)
                    or checker.check_reaction("Stille reaction_vinyl", rsmi)
                    or checker.check_reaction("Stille reaction_benzyl", rsmi)
                    or checker.check_reaction("Stille reaction_allyl", rsmi)
                    or checker.check_reaction("Hiyama-Denmark Coupling", rsmi)
                    or checker.check_reaction("Kumada cross-coupling", rsmi)
                    or checker.check_reaction("Aryllithium cross-coupling", rsmi)
                ):

                    print(f"Detected late-stage cross-coupling at depth {depth}")
                    has_late_stage_coupling = True

                # If no specific reaction match, check for characteristic functional groups
                if not has_late_stage_coupling:
                    reactants = rsmi.split(">")[0].split(".")

                    # Check for boronic compounds
                    has_boronic = any(
                        checker.check_fg("Boronic acid", r) or checker.check_fg("Boronic ester", r)
                        for r in reactants
                    )

                    # Check for halides
                    has_halide = any(
                        checker.check_fg("Aromatic halide", r)
                        or checker.check_fg("Primary halide", r)
                        or checker.check_fg("Secondary halide", r)
                        or checker.check_fg("Tertiary halide", r)
                        or checker.check_fg("Alkenyl halide", r)
                        for r in reactants
                    )

                    # Check for other organometallic compounds
                    has_organometallic = any(
                        checker.check_fg("Magnesium halide", r) or "Sn" in r or "Zn" in r
                        for r in reactants
                    )

                    if (has_boronic or has_organometallic) and has_halide:
                        print(
                            f"Detected potential late-stage cross-coupling at depth {depth} based on functional groups"
                        )
                        has_late_stage_coupling = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_late_stage_coupling

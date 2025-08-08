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
    Detects if the synthesis route employs a Suzuki coupling strategy,
    looking for boronic acid derivatives and halogenated aromatics as reactants.
    """
    suzuki_coupling_found = False

    def dfs_traverse(node):
        nonlocal suzuki_coupling_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a Suzuki coupling reaction using the checker function
            if (
                checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with sulfonic esters", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                or checker.check_reaction("Suzuki", rsmi)
            ):

                # Extract depth information safely
                depth = 999  # Default depth if not found
                if "ID" in node["metadata"]:
                    depth_match = re.search(r"Depth: (\d+)", node["metadata"]["ID"])
                    if depth_match:
                        depth = int(depth_match.group(1))

                print(f"Suzuki coupling detected at depth {depth}")
                suzuki_coupling_found = True

            # Fallback method if checker doesn't identify the reaction
            if not suzuki_coupling_found:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                has_boronic_acid = False
                has_boronic_ester = False
                has_halogen = False

                for r in reactants:
                    if checker.check_fg("Boronic acid", r):
                        has_boronic_acid = True
                    if checker.check_fg("Boronic ester", r):
                        has_boronic_ester = True
                    if checker.check_fg("Aromatic halide", r):
                        has_halogen = True

                if (has_boronic_acid or has_boronic_ester) and has_halogen:
                    # Extract depth information safely
                    depth = 999  # Default depth if not found
                    if "ID" in node["metadata"]:
                        depth_match = re.search(r"Depth: (\d+)", node["metadata"]["ID"])
                        if depth_match:
                            depth = int(depth_match.group(1))

                    print(f"Suzuki coupling detected (fallback method) at depth {depth}")
                    suzuki_coupling_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_coupling_found

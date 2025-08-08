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
    This function detects the use of a multi-component reaction (like Ugi reaction)
    in the late stage of the synthesis.
    """
    late_stage_mcr = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_mcr

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if this is a multi-component reaction (3+ reactants)
                if len(reactants) >= 3:
                    # Check for known MCR reaction types
                    is_ugi = checker.check_reaction("Ugi reaction", rsmi)
                    is_a3_coupling = checker.check_reaction("A3 coupling", rsmi)

                    # Check for isocyanide component (characteristic of Ugi reaction)
                    has_isocyanide = any(checker.check_fg("Isocyanide", r) for r in reactants)

                    if is_ugi or is_a3_coupling or has_isocyanide:
                        late_stage_mcr = True
                        print(f"Detected late-stage multi-component reaction at depth {depth}")
                        print(f"Reaction SMILES: {rsmi}")
                        print(f"Number of reactants: {len(reactants)}")
                        if is_ugi:
                            print("Identified as Ugi reaction")
                        elif is_a3_coupling:
                            print("Identified as A3 coupling")
                        elif has_isocyanide:
                            print("Contains isocyanide component (likely Ugi-type reaction)")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return late_stage_mcr

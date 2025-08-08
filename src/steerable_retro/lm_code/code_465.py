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
    This function detects a synthetic strategy involving late-stage N-alkylation
    to join complex fragments.
    """
    has_late_stage_n_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_n_alkylation

        # Focus on reactions at depths 0-1 (late stage)
        if node["type"] == "reaction" and depth <= 1:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an N-alkylation reaction
                is_n_alkylation = checker.check_reaction(
                    "N-alkylation of primary amines with alkyl halides", rsmi
                ) or checker.check_reaction(
                    "N-alkylation of secondary amines with alkyl halides", rsmi
                )

                if not is_n_alkylation:
                    # Try checking for specific functional groups if reaction check fails
                    has_primary_amine = any(checker.check_fg("Primary amine", r) for r in reactants)
                    has_secondary_amine = any(
                        checker.check_fg("Secondary amine", r) for r in reactants
                    )
                    has_alkyl_halide = any(
                        checker.check_fg("Primary halide", r)
                        or checker.check_fg("Secondary halide", r)
                        or checker.check_fg("Tertiary halide", r)
                        for r in reactants
                    )
                    has_tertiary_amine = checker.check_fg("Tertiary amine", product)

                    is_n_alkylation = (
                        (has_primary_amine or has_secondary_amine)
                        and has_alkyl_halide
                        and has_tertiary_amine
                    )

                if is_n_alkylation:
                    # Check if we have at least 2 reactants (fragment joining)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    if len(reactant_mols) >= 2:
                        # Count heavy atoms to ensure these are complex fragments
                        reactant_heavy_atoms = [
                            sum(1 for atom in r.GetAtoms() if atom.GetAtomicNum() > 1)
                            for r in reactant_mols
                            if r
                        ]
                        if any(count >= 10 for count in reactant_heavy_atoms):
                            print(
                                f"Detected late-stage N-alkylation joining complex fragments at depth {depth}"
                            )
                            has_late_stage_n_alkylation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return has_late_stage_n_alkylation

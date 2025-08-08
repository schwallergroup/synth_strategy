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
    This function detects if the synthesis route involves a benzyl halide intermediate.
    A benzyl halide intermediate is a compound that contains both a benzene ring and a primary halide,
    and is used as a reactant in a reaction step.
    """
    benzyl_halide_nodes = []
    benzyl_halide_as_reactant = False

    def dfs_traverse(node, depth=0):
        nonlocal benzyl_halide_as_reactant

        if node["type"] == "mol":
            # Check if the molecule is a benzyl halide (has benzene ring and primary halide)
            if checker.check_ring("benzene", node["smiles"]) and checker.check_fg(
                "Primary halide", node["smiles"]
            ):
                print(f"Found benzyl halide at depth {depth}: {node['smiles']}")
                benzyl_halide_nodes.append((node["smiles"], depth))

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Check if any benzyl halide is used as a reactant in this reaction
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            for reactant in reactants:
                if checker.check_ring("benzene", reactant) and checker.check_fg(
                    "Primary halide", reactant
                ):
                    print(f"Benzyl halide used as reactant in reaction: {rsmi}")
                    benzyl_halide_as_reactant = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we found benzyl halides at different depths AND at least one is used as a reactant
    has_benzyl_halide_intermediate = (
        len(set([depth for _, depth in benzyl_halide_nodes])) >= 1 and benzyl_halide_as_reactant
    )

    print(f"Benzyl halide intermediate strategy detected: {has_benzyl_halide_intermediate}")
    print(f"Benzyl halide nodes found: {benzyl_halide_nodes}")
    return has_benzyl_halide_intermediate

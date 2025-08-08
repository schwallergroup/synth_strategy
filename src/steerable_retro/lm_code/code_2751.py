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
    Detects a linear synthesis strategy that includes a Suzuki coupling step.
    """
    found_suzuki = False
    is_linear = True
    main_path_nodes = []

    def dfs_traverse(node, depth=0, path=None):
        nonlocal found_suzuki, is_linear, main_path_nodes

        if path is None:
            path = []

        # Track this node in the current path
        current_path = path + [node]

        # If this is a leaf node (starting material or no children)
        if node["type"] == "mol" and (
            node.get("in_stock", False) or len(node.get("children", [])) == 0
        ):
            # If this is the longest path so far, update main_path_nodes
            if len(current_path) > len(main_path_nodes):
                main_path_nodes = current_path

        # Check if this node branches (non-linear)
        if node["type"] == "mol" and len(node.get("children", [])) > 1:
            # We'll allow some branching, but track it
            print(
                f"Potential branch found at depth {depth} with {len(node.get('children', []))} children"
            )

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check for Suzuki coupling using the checker function
            suzuki_types = [
                "Suzuki coupling with boronic acids",
                "Suzuki coupling with boronic esters",
                "Suzuki coupling with boronic acids OTf",
                "Suzuki coupling with boronic esters OTf",
                "Suzuki coupling with sulfonic esters",
                "Suzuki",  # Generic Suzuki check
            ]

            for suzuki_type in suzuki_types:
                if checker.check_reaction(suzuki_type, rsmi):
                    found_suzuki = True
                    print(f"Found {suzuki_type} at depth {depth}")
                    break

            # If we haven't found a Suzuki reaction yet, check for characteristic functional groups
            if not found_suzuki:
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for boronic acids/esters in reactants
                    has_boronic = False
                    has_aryl_halide = False

                    for reactant in reactants:
                        if checker.check_fg("Boronic acid", reactant) or checker.check_fg(
                            "Boronic ester", reactant
                        ):
                            has_boronic = True
                            print(f"Found boronic acid/ester in reactant: {reactant}")

                        if (
                            checker.check_fg("Aromatic halide", reactant)
                            or checker.check_fg("Aromatic chloride", reactant)
                            or checker.check_fg("Aromatic bromide", reactant)
                            or checker.check_fg("Aromatic iodide", reactant)
                            or checker.check_fg("Triflate", reactant)
                        ):
                            has_aryl_halide = True
                            print(f"Found aryl halide/triflate in reactant: {reactant}")

                    # If we have both boronic compound and aryl halide, it's likely a Suzuki coupling
                    if has_boronic and has_aryl_halide:
                        found_suzuki = True
                        print(f"Identified Suzuki coupling by functional groups at depth {depth}")
                except Exception as e:
                    print(f"Error analyzing reaction components: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

    # Start traversal
    dfs_traverse(route)

    # Check if the main synthetic pathway is linear
    # We'll define "linear" as having no more than one reaction child per molecule node
    # in the main synthetic pathway
    if main_path_nodes:
        for i, node in enumerate(main_path_nodes):
            if node["type"] == "mol" and i < len(main_path_nodes) - 1:
                # Count reaction children that are in the main path
                reaction_children_in_main_path = 0
                for child in node.get("children", []):
                    if child in main_path_nodes:
                        reaction_children_in_main_path += 1

                if reaction_children_in_main_path > 1:
                    is_linear = False
                    print(
                        f"Non-linear main pathway: molecule has {reaction_children_in_main_path} reaction children in main path"
                    )

    strategy_present = found_suzuki and is_linear
    print(f"Linear synthesis with Suzuki coupling strategy: {strategy_present}")
    print(f"- Found Suzuki coupling: {found_suzuki}")
    print(f"- Is linear synthesis: {is_linear}")
    print(f"- Main path length: {len(main_path_nodes)}")

    return strategy_present

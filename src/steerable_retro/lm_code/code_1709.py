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
    Detects if the synthetic route involves multiple cross-coupling reactions
    (likely Suzuki couplings) for building the molecular scaffold.
    """
    cross_coupling_count = 0
    nodes_visited = 0

    def dfs_traverse(node, depth=0):
        nonlocal cross_coupling_count, nodes_visited
        nodes_visited += 1

        print(f"Visiting node at depth {depth}: {node.get('type')} - {node.get('smiles', '')}")

        if node["type"] == "reaction":
            if node.get("metadata", {}).get("rsmi"):
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction SMILES: {rsmi}")

                # Extract reactants and product
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for various cross-coupling reactions using the checker
                cross_coupling_reactions = [
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic esters",
                    "Suzuki coupling with boronic acids OTf",
                    "Suzuki coupling with boronic esters OTf",
                    "Suzuki coupling with sulfonic esters",
                    "Negishi coupling",
                    "Stille reaction_vinyl",
                    "Stille reaction_aryl",
                    "Stille reaction_benzyl",
                    "Stille reaction_allyl",
                    "Stille reaction_vinyl OTf",
                    "Stille reaction_aryl OTf",
                    "Stille reaction_benzyl OTf",
                    "Stille reaction_allyl OTf",
                    "Stille reaction_other",
                    "Stille reaction_other OTf",
                    "Heck terminal vinyl",
                    "Oxidative Heck reaction",
                    "Oxidative Heck reaction with vinyl ester",
                    "Heck reaction with vinyl ester and amine",
                    "Sonogashira acetylene_aryl halide",
                    "Sonogashira alkyne_aryl halide",
                    "Sonogashira acetylene_aryl OTf",
                    "Sonogashira alkyne_aryl OTf",
                    "Sonogashira acetylene_alkenyl halide",
                    "Sonogashira alkyne_alkenyl halide",
                    "Sonogashira acetylene_alkenyl OTf",
                    "Sonogashira alkyne_alkenyl OTf",
                    "Sonogashira acetylene_acyl halide",
                    "Sonogashira alkyne_acyl halide",
                    "Kumada cross-coupling",
                    "Hiyama-Denmark Coupling",
                    "Aryllithium cross-coupling",
                    "Buchwald-Hartwig",
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "{Suzuki}",  # Additional pattern from the list
                ]

                # First try to match using reaction checkers
                is_cross_coupling = False
                for reaction_type in cross_coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        cross_coupling_count += 1
                        print(f"Detected cross-coupling reaction: {reaction_type} in {rsmi}")
                        is_cross_coupling = True
                        break  # Count each reaction only once

                # If no match, check for characteristic functional groups and catalysts
                if not is_cross_coupling:
                    # Check for palladium catalyst in the reagents section
                    reagents = rsmi.split(">")[1].split(".")
                    has_pd_catalyst = any("Pd" in reagent for reagent in reagents)

                    # Check for boronic acid/ester in reactants
                    has_boronic = any(
                        checker.check_fg("Boronic acid", reactant)
                        or checker.check_fg("Boronic ester", reactant)
                        for reactant in reactants
                    )

                    # Check for aryl halide in reactants
                    has_aryl_halide = any(
                        checker.check_fg("Aromatic halide", reactant) for reactant in reactants
                    )

                    # Check for other metal-containing reactants
                    has_metal_reactant = any(
                        checker.check_fg("Magnesium halide", reactant)
                        or "Sn" in reactant
                        or "Zn" in reactant
                        for reactant in reactants
                    )

                    # If we have characteristic features of cross-coupling
                    if has_pd_catalyst and (has_boronic or has_metal_reactant) and has_aryl_halide:
                        cross_coupling_count += 1
                        print(f"Detected cross-coupling reaction by functional groups in {rsmi}")
                        print(f"  Pd catalyst: {has_pd_catalyst}")
                        print(f"  Boronic acid/ester: {has_boronic}")
                        print(f"  Aryl halide: {has_aryl_halide}")
                        print(f"  Metal reactant: {has_metal_reactant}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Total nodes visited: {nodes_visited}")
    print(f"Total cross-coupling reactions found: {cross_coupling_count}")
    return cross_coupling_count >= 2

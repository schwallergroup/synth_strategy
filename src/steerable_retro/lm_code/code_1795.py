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
    Detects a synthesis that utilizes halogens (Cl, Br, I) as synthetic handles
    throughout multiple steps of the route.
    """
    # Track halogen presence and transformations at different depths
    halogen_depths = set()
    halogen_reaction_depths = set()

    def dfs_traverse(node, depth=0):
        # For molecule nodes, check for halogen presence
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for various halogen functional groups
            halogen_present = (
                checker.check_fg("Aromatic halide", mol_smiles)
                or checker.check_fg("Primary halide", mol_smiles)
                or checker.check_fg("Secondary halide", mol_smiles)
                or checker.check_fg("Tertiary halide", mol_smiles)
                or checker.check_fg("Alkenyl halide", mol_smiles)
                or checker.check_fg("Haloalkyne", mol_smiles)
            )

            if halogen_present:
                halogen_depths.add(depth)
                print(f"Found halogen in molecule at depth {depth}: {mol_smiles}")

        # For reaction nodes, check if the reaction involves halogen transformation
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for reactions involving halogens
            halogen_reaction = (
                checker.check_reaction("Aromatic fluorination", rxn_smiles)
                or checker.check_reaction("Aromatic chlorination", rxn_smiles)
                or checker.check_reaction("Aromatic bromination", rxn_smiles)
                or checker.check_reaction("Aromatic iodination", rxn_smiles)
                or checker.check_reaction("Chlorination", rxn_smiles)
                or checker.check_reaction("Fluorination", rxn_smiles)
                or checker.check_reaction("Iodination", rxn_smiles)
                or checker.check_reaction("Bromination", rxn_smiles)
                or checker.check_reaction(
                    "Aromatic substitution of bromine by chlorine", rxn_smiles
                )
                or checker.check_reaction("Aromatic dehalogenation", rxn_smiles)
                or checker.check_reaction("Dehalogenation", rxn_smiles)
                or checker.check_reaction("Finkelstein reaction", rxn_smiles)
                or checker.check_reaction("Halodeboronation of boronic acids", rxn_smiles)
                or checker.check_reaction("Halodeboronation of boronic esters", rxn_smiles)
            )

            if halogen_reaction:
                halogen_reaction_depths.add(depth)
                print(f"Found halogen-involving reaction at depth {depth}: {rxn_smiles}")

            # Check for reactions that commonly use halogens as leaving groups
            nucleophilic_reactions = (
                checker.check_reaction("Suzuki coupling with boronic acids", rxn_smiles)
                or checker.check_reaction("Suzuki coupling with boronic esters", rxn_smiles)
                or checker.check_reaction("Negishi coupling", rxn_smiles)
                or checker.check_reaction("Heck terminal vinyl", rxn_smiles)
                or checker.check_reaction("Sonogashira acetylene_aryl halide", rxn_smiles)
                or checker.check_reaction("Sonogashira alkyne_aryl halide", rxn_smiles)
                or checker.check_reaction(
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rxn_smiles
                )
                or checker.check_reaction(
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rxn_smiles
                )
            )

            if nucleophilic_reactions:
                # Check if reactants contain halogens
                try:
                    reactants = rxn_smiles.split(">")[0].split(".")
                    for reactant in reactants:
                        if (
                            checker.check_fg("Aromatic halide", reactant)
                            or checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                        ):
                            halogen_reaction_depths.add(depth)
                            print(
                                f"Found coupling reaction using halogen at depth {depth}: {rxn_smiles}"
                            )
                            break
                except Exception as e:
                    print(f"Error checking reactants: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if halogens are used in at least 2 different depths
    # Consider both molecule presence and reaction involvement
    all_halogen_depths = halogen_depths.union(halogen_reaction_depths)
    result = len(all_halogen_depths) >= 2

    print(f"Halogen depths: {halogen_depths}")
    print(f"Halogen reaction depths: {halogen_reaction_depths}")
    print(f"Halogen-mediated synthesis across multiple steps: {result}")

    return result

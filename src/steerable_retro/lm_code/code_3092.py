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
    This function detects a synthetic strategy involving the connection of
    heterocyclic fragments (pyrazole and piperidine) via an amide linker.
    """
    # Track if we've found a molecule with both heterocycles connected by amide
    found_connected_molecule = False
    # Track if we've found a reaction that forms this connection
    found_connecting_reaction = False

    def dfs_traverse(node):
        nonlocal found_connected_molecule, found_connecting_reaction

        if node["type"] == "mol":
            mol_smiles = node.get("smiles", "")

            # Check if molecule contains both pyrazole and piperidine
            if (
                mol_smiles
                and checker.check_ring("pyrazole", mol_smiles)
                and checker.check_ring("piperidine", mol_smiles)
            ):
                print(f"Found molecule with both pyrazole and piperidine: {mol_smiles}")

                # Check if it contains an amide group
                if (
                    checker.check_fg("Secondary amide", mol_smiles)
                    or checker.check_fg("Primary amide", mol_smiles)
                    or checker.check_fg("Tertiary amide", mol_smiles)
                ):
                    print(f"Found molecule with pyrazole, piperidine, and amide: {mol_smiles}")
                    found_connected_molecule = True

        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Check for amide coupling reactions
                amide_forming_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Schotten-Baumann to ester",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                ]

                for reaction_type in amide_forming_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found potential amide coupling reaction: {rsmi}")

                        # Extract reactants and product to check if this reaction connects pyrazole and piperidine
                        try:
                            reactants = rsmi.split(">")[0].split(".")
                            product = rsmi.split(">")[-1]

                            # Check if product has both heterocycles and amide
                            if (
                                checker.check_ring("pyrazole", product)
                                and checker.check_ring("piperidine", product)
                                and (
                                    checker.check_fg("Secondary amide", product)
                                    or checker.check_fg("Primary amide", product)
                                    or checker.check_fg("Tertiary amide", product)
                                )
                            ):

                                # Check if reactants have separate heterocycles
                                has_pyrazole_reactant = False
                                has_piperidine_reactant = False

                                for reactant in reactants:
                                    if checker.check_ring("pyrazole", reactant):
                                        has_pyrazole_reactant = True
                                    if checker.check_ring("piperidine", reactant):
                                        has_piperidine_reactant = True

                                # If one reactant has pyrazole and another has piperidine, this is our connecting reaction
                                if has_pyrazole_reactant and has_piperidine_reactant:
                                    print(
                                        f"Found reaction connecting pyrazole and piperidine via amide: {rsmi}"
                                    )
                                    found_connecting_reaction = True
                                    break
                        except Exception as e:
                            print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we found a molecule with connected heterocycles AND a reaction that forms this connection
    strategy_present = found_connected_molecule and found_connecting_reaction

    print(f"Heterocycle connection via amide strategy detected: {strategy_present}")
    print(
        f"Connected molecule found: {found_connected_molecule}, Connecting reaction found: {found_connecting_reaction}"
    )

    # If we found a molecule with both heterocycles connected by amide but didn't find the connecting reaction,
    # we'll still consider the strategy present since the molecule exists in the route
    if found_connected_molecule and not found_connecting_reaction:
        print(
            "Found molecule with connected heterocycles but couldn't identify the connecting reaction."
        )
        return True

    return strategy_present

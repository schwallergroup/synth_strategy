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
    This function detects if the route involves late-stage incorporation of a morpholine group.
    """
    morpholine_added = False
    late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal morpholine_added, late_stage

        if node["type"] == "reaction" and depth <= 2:  # Late stage (expanded definition)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check if product contains morpholine
                if checker.check_ring("morpholine", product):
                    print(f"Product contains morpholine at depth {depth}")

                    # Check if morpholine is being incorporated (present in reactants as separate entity)
                    morpholine_in_reactants = False
                    non_morpholine_in_reactants = False

                    for reactant in reactants:
                        if checker.check_ring("morpholine", reactant):
                            # If the reactant is primarily just morpholine (or a simple derivative)
                            if len(Chem.MolFromSmiles(reactant).GetAtoms()) <= 15:
                                morpholine_in_reactants = True
                                print(f"Found morpholine-containing reactant: {reactant}")
                        else:
                            non_morpholine_in_reactants = True

                    # Check for relevant reaction types that would incorporate morpholine
                    is_incorporation_reaction = False

                    # Check common reaction types for morpholine incorporation
                    if checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                    ):
                        is_incorporation_reaction = True
                        print(f"Detected morpholine incorporation reaction: Sulfonamide synthesis")
                    elif checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    ) or checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    ):
                        is_incorporation_reaction = True
                        print(f"Detected morpholine incorporation reaction: N-alkylation")
                    elif checker.check_reaction("Williamson Ether Synthesis", rsmi):
                        is_incorporation_reaction = True
                        print(
                            f"Detected morpholine incorporation reaction: Williamson Ether Synthesis"
                        )
                    elif checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                    ) or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                    ):
                        is_incorporation_reaction = True
                        print(f"Detected morpholine incorporation reaction: N-arylation")
                    elif checker.check_reaction("sulfon_amide", rsmi):
                        is_incorporation_reaction = True
                        print(f"Detected morpholine incorporation reaction: Sulfonamide formation")

                    # Check for sulfonylation specifically
                    if morpholine_in_reactants and any(
                        "S(=O)(=O)Cl" in r or "S(Cl)(=O)=O" in r for r in reactants
                    ):
                        is_incorporation_reaction = True
                        print(f"Detected sulfonylation of morpholine at depth {depth}")

                    # If we have morpholine in reactants, non-morpholine in reactants, and it's an incorporation reaction
                    if morpholine_in_reactants and non_morpholine_in_reactants:
                        morpholine_added = True
                        late_stage = True
                        print(f"Late-stage morpholine incorporation confirmed at depth {depth}")

                    # Alternative detection: if product has morpholine but no reactant has it
                    elif not morpholine_in_reactants and non_morpholine_in_reactants:
                        # This means morpholine was formed in the reaction
                        morpholine_added = True
                        late_stage = True
                        print(f"Late-stage morpholine formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: morpholine_added={morpholine_added}, late_stage={late_stage}")
    return morpholine_added and late_stage

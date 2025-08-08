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
    Detects if the synthesis follows a linear strategy (no convergent steps)
    where each reaction has only one non-reagent reactant.
    """
    is_linear = True

    # List of common reagent functional groups
    reagent_fg = [
        "Triflate",
        "Tosylate",
        "Mesylate",
        "Boronic acid",
        "Boronic ester",
        "Magnesium halide",
        "Zinc halide",
        "Tin",
        "Alkyl lithium",
        "Silane",
        "Acyl halide",
        "Anhydride",
    ]

    # Common protecting groups and their reagents
    protecting_reagents = [
        "Boc",
        "TMS ether protective group",
        "Silyl protective group",
        "Acetal/Ketal",
    ]

    # Reactions that typically have special handling
    coupling_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Sonogashira alkyne_aryl halide",
        "Sonogashira acetylene_aryl halide",
        "Buchwald-Hartwig",
        "Negishi coupling",
        "Stille reaction",
        "Heck terminal vinyl",
        "Ullmann-Goldberg Substitution",
    ]

    protection_reactions = [
        "Boc amine protection",
        "Alcohol protection with silyl ethers",
        "Protection of carboxylic acid",
        "Aldehyde or ketone acetalization",
    ]

    deprotection_reactions = [
        "Boc amine deprotection",
        "Alcohol deprotection from silyl ethers",
        "Deprotection of carboxylic acid",
        "Acetal hydrolysis",
    ]

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                print(f"Processing reaction: {rsmi}")

                # Check if this is a coupling, protection, or deprotection reaction
                is_coupling = any(checker.check_reaction(rxn, rsmi) for rxn in coupling_reactions)
                is_protection = any(
                    checker.check_reaction(rxn, rsmi) for rxn in protection_reactions
                )
                is_deprotection = any(
                    checker.check_reaction(rxn, rsmi) for rxn in deprotection_reactions
                )

                print(
                    f"Reaction classification - Coupling: {is_coupling}, Protection: {is_protection}, Deprotection: {is_deprotection}"
                )

                # Count substantial reactants (excluding small reagents)
                substantial_reactants = []
                for r in reactants:
                    # Skip empty strings
                    if not r:
                        continue

                    try:
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            # Count atoms excluding H
                            atom_count = mol.GetNumHeavyAtoms()

                            # Check if this is likely a reagent
                            is_reagent = False

                            # Small molecules are typically reagents
                            if atom_count <= 3:
                                is_reagent = True
                                print(
                                    f"Identified small molecule reagent: {r} with {atom_count} heavy atoms"
                                )

                            # Check for common reagent functional groups
                            for fg in reagent_fg:
                                if checker.check_fg(fg, r):
                                    is_reagent = True
                                    print(f"Identified reagent with {fg} functional group: {r}")
                                    break

                            # Check for protecting group reagents
                            for pg in protecting_reagents:
                                if checker.check_fg(pg, r):
                                    is_reagent = True
                                    print(f"Identified protecting group reagent: {r}")
                                    break

                            # Special handling for coupling reactions
                            if is_coupling and not is_reagent and atom_count > 3:
                                # For coupling reactions, check if this is a coupling partner
                                coupling_groups = [
                                    "Boronic acid",
                                    "Boronic ester",
                                    "Magnesium halide",
                                    "Zinc halide",
                                    "Tin",
                                    "Alkyne",
                                    "Aromatic halide",
                                ]

                                # If we already have one substantial reactant, this one might be a coupling partner
                                if len(substantial_reactants) > 0:
                                    for fg in coupling_groups:
                                        if checker.check_fg(fg, r):
                                            is_reagent = True
                                            print(f"Identified coupling partner with {fg}: {r}")
                                            break

                            # Special handling for protection reactions
                            if is_protection and not is_reagent and atom_count > 3:
                                # Check if this is a protecting group reagent
                                if any(checker.check_fg(pg, r) for pg in protecting_reagents):
                                    is_reagent = True
                                    print(f"Identified protection reagent: {r}")

                            # Special handling for Boc anhydride and similar large reagents
                            if not is_reagent and "OC(=O)OC(=O)" in r:
                                is_reagent = True
                                print(f"Identified anhydride reagent: {r}")

                            # Consider molecules with more than 3 heavy atoms as substantial if not a reagent
                            if not is_reagent:
                                substantial_reactants.append(r)
                                print(
                                    f"Identified substantial reactant: {r} with {atom_count} heavy atoms"
                                )
                    except Exception as e:
                        print(f"Error processing reactant {r}: {e}")

                # If more than one substantial reactant, it's not a linear synthesis
                if len(substantial_reactants) > 1:
                    is_linear = False
                    print(
                        f"Found convergent step with {len(substantial_reactants)} substantial reactants: {substantial_reactants}"
                    )
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Linear synthesis detection result: {is_linear}")
    return is_linear

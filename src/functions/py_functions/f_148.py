#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    Detects if the synthesis route follows a protection-coupling-deprotection pattern.
    """
    # Track protection, coupling, and deprotection events in order
    events = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction at depth {depth}: {rsmi}")
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for protection reactions
                protection_reactions = [
                    "Boc amine protection",
                    "Boc amine protection explicit",
                    "Boc amine protection with Boc anhydride",
                    "Boc amine protection (ethyl Boc)",
                    "Boc amine protection of secondary amine",
                    "Boc amine protection of primary amine",
                    "Alcohol protection with silyl ethers",
                    "Protection of carboxylic acid",
                ]

                # Check for coupling reactions
                coupling_reactions = [
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic esters",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "Sonogashira acetylene_aryl halide",
                    "Sonogashira alkyne_aryl halide",
                    "Heck terminal vinyl",
                    "Negishi coupling",
                    "Stille reaction_aryl",
                    "Ullmann condensation",
                    "Goldberg coupling",
                    "{Suzuki}",
                    "{N-arylation_heterocycles}",
                    "{Buchwald-Hartwig}",
                ]

                # Check for deprotection reactions
                deprotection_reactions = [
                    "Boc amine deprotection",
                    "Boc amine deprotection of guanidine",
                    "Boc amine deprotection to NH-NH2",
                    "Alcohol deprotection from silyl ethers",
                    "Alcohol deprotection from silyl ethers (double)",
                    "Alcohol deprotection from silyl ethers (diol)",
                    "Deprotection of carboxylic acid",
                    "Ester saponification (methyl deprotection)",
                    "Ester saponification (alkyl deprotection)",
                    "Hydroxyl benzyl deprotection",
                    "Carboxyl benzyl deprotection",
                    "COOH ethyl deprotection",
                    "Phthalimide deprotection",
                    "Tert-butyl deprotection of amine",
                    "TMS deprotection from alkyne",
                    "N-glutarimide deprotection",
                ]

                # Check for protection event
                protection_found = False
                for protection_rxn in protection_reactions:
                    if checker.check_reaction(protection_rxn, rsmi):
                        events.append(("protection", depth, rsmi))
                        print(f"Found protection at depth {depth}: {protection_rxn}")
                        protection_found = True
                        break

                # Fallback detection for protection using functional groups
                if not protection_found:
                    # Check for Boc protection
                    if any(
                        checker.check_fg("Primary amine", r) for r in reactants
                    ) and checker.check_fg("Boc", product):
                        events.append(("protection", depth, rsmi))
                        print(f"Found Boc protection via FG analysis at depth {depth}")
                    # Check for silyl protection
                    elif any(
                        checker.check_fg("Primary alcohol", r) for r in reactants
                    ) and checker.check_fg("Silyl protective group", product):
                        events.append(("protection", depth, rsmi))
                        print(
                            f"Found silyl protection via FG analysis at depth {depth}"
                        )

                # Check for coupling event
                coupling_found = False
                for coupling_rxn in coupling_reactions:
                    if checker.check_reaction(coupling_rxn, rsmi):
                        events.append(("coupling", depth, rsmi))
                        print(f"Found coupling at depth {depth}: {coupling_rxn}")
                        coupling_found = True
                        break

                # Fallback detection for coupling using functional groups
                if not coupling_found:
                    # Check for C-N bond formation (Buchwald-Hartwig/N-arylation)
                    if any(
                        checker.check_fg("Aromatic halide", r) for r in reactants
                    ) and any(
                        checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Primary amine", r)
                        or checker.check_fg("Tertiary amine", r)
                        for r in reactants
                    ):
                        events.append(("coupling", depth, rsmi))
                        print(f"Found C-N coupling via FG analysis at depth {depth}")
                    # Check for C-C bond formation (common in couplings)
                    elif any(
                        checker.check_fg("Aromatic halide", r) for r in reactants
                    ) and any(
                        checker.check_fg("Boronic acid", r)
                        or checker.check_fg("Boronic ester", r)
                        for r in reactants
                    ):
                        events.append(("coupling", depth, rsmi))
                        print(f"Found C-C coupling via FG analysis at depth {depth}")

                # Check for deprotection event
                deprotection_found = False
                for deprotection_rxn in deprotection_reactions:
                    if checker.check_reaction(deprotection_rxn, rsmi):
                        events.append(("deprotection", depth, rsmi))
                        print(
                            f"Found deprotection at depth {depth}: {deprotection_rxn}"
                        )
                        deprotection_found = True
                        break

                # Fallback detection for deprotection using functional groups
                if not deprotection_found:
                    # Check for Boc deprotection
                    if any(checker.check_fg("Boc", r) for r in reactants) and (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                    ):
                        events.append(("deprotection", depth, rsmi))
                        print(
                            f"Found Boc deprotection via FG analysis at depth {depth}"
                        )
                    # Check for silyl deprotection
                    elif any(
                        checker.check_fg("Silyl protective group", r) for r in reactants
                    ) and checker.check_fg("Primary alcohol", product):
                        events.append(("deprotection", depth, rsmi))
                        print(
                            f"Found silyl deprotection via FG analysis at depth {depth}"
                        )
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check for protection-coupling-deprotection pattern
    # Sort events by depth to ensure chronological order in the synthesis direction
    # Note: In retrosynthesis, higher depth = earlier in actual synthesis
    events.sort(key=lambda x: x[1], reverse=True)

    print(f"Sorted events: {events}")

    if not events:
        print("No protection, coupling, or deprotection events found")
        return False

    # Extract depths for each event type
    protection_depths = [d for t, d, _ in events if t == "protection"]
    coupling_depths = [d for t, d, _ in events if t == "coupling"]
    deprotection_depths = [d for t, d, _ in events if t == "deprotection"]

    # Check if we have all three event types
    if not (protection_depths and coupling_depths and deprotection_depths):
        print(
            f"Missing event types. Protection: {bool(protection_depths)}, Coupling: {bool(coupling_depths)}, Deprotection: {bool(deprotection_depths)}"
        )
        return False

    # Check if they occur in the correct order in the synthesis direction
    # In synthesis direction (reverse of retrosynthesis), we want:
    # protection (highest depth) -> coupling -> deprotection (lowest depth)
    if max(protection_depths) > max(coupling_depths) > max(deprotection_depths):
        print("Found protection-coupling-deprotection pattern in correct order")
        return True

    # Alternative check: look for the pattern in the sequence of events
    found_protection = False
    found_coupling = False
    found_deprotection = False

    for event_type, depth, _ in events:
        if event_type == "protection" and not found_protection:
            found_protection = True
        elif event_type == "coupling" and found_protection and not found_coupling:
            found_coupling = True
        elif event_type == "deprotection" and found_coupling and not found_deprotection:
            found_deprotection = True

    if found_protection and found_coupling and found_deprotection:
        print("Found protection-coupling-deprotection pattern in sequential order")
        return True

    return False

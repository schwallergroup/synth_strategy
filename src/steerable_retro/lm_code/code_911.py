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
    This function detects if Boc protection is maintained through multiple steps.
    A proper Boc protection strategy involves:
    1. Adding Boc protection early in the synthesis
    2. Maintaining the Boc group through intermediate steps
    3. Removing the Boc group in a late stage

    Note: In retrosynthetic analysis, higher depth = earlier stage, lower depth = later stage
    """
    # Track Boc-protected molecules and protection/deprotection reactions
    protected_molecules = []
    protection_reactions = []
    deprotection_reactions = []

    # Track depth of each reaction for determining early vs late stage
    reaction_depths = {}
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            # Check if molecule has Boc group
            if checker.check_fg("Boc", node["smiles"]):
                protected_molecules.append((depth, node["smiles"]))
                print(f"Found Boc-protected molecule at depth {depth}: {node['smiles'][:30]}...")

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Store reaction depth for later analysis
                reaction_depths[rsmi] = depth

                # Check if this is a Boc protection reaction
                if (
                    checker.check_reaction("Boc amine protection", rsmi)
                    or checker.check_reaction("Boc amine protection explicit", rsmi)
                    or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                    or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                    or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                    or checker.check_reaction("Boc amine protection of primary amine", rsmi)
                ):
                    protection_reactions.append((depth, rsmi))
                    print(f"Found Boc protection reaction at depth {depth}: {rsmi[:50]}...")

                # Check if this is a Boc deprotection reaction
                if (
                    checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                    or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                ):
                    deprotection_reactions.append((depth, rsmi))
                    print(f"Found Boc deprotection reaction at depth {depth}: {rsmi[:50]}...")

                # Alternative check: see if Boc group appears or disappears
                has_boc_in_reactants = any(checker.check_fg("Boc", r) for r in reactants)
                has_boc_in_product = checker.check_fg("Boc", product)

                # In retrosynthetic analysis, if product has Boc but reactants don't,
                # it's a deprotection reaction when going forward
                if has_boc_in_reactants and not has_boc_in_product:
                    if not any(rsmi == r[1] for r in deprotection_reactions):
                        deprotection_reactions.append((depth, rsmi))
                        print(f"Found implicit Boc deprotection at depth {depth}: {rsmi[:50]}...")

                # In retrosynthetic analysis, if reactants have Boc but product doesn't,
                # it's a protection reaction when going forward
                if has_boc_in_product and not has_boc_in_reactants:
                    if not any(rsmi == r[1] for r in protection_reactions):
                        protection_reactions.append((depth, rsmi))
                        print(f"Found implicit Boc protection at depth {depth}: {rsmi[:50]}...")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Analyze the protection strategy
    print(f"Found {len(protection_reactions)} Boc protection reactions")
    print(f"Found {len(deprotection_reactions)} Boc deprotection reactions")
    print(f"Found {len(protected_molecules)} Boc-protected molecules")

    # No protection or deprotection found
    if not protection_reactions and not deprotection_reactions:
        print("No Boc protection/deprotection reactions found")
        return False

    # If we found deprotection but no protection, try to infer protection
    if deprotection_reactions and not protection_reactions:
        print("Found deprotection but no explicit protection, checking for protected molecules")
        # If we have protected molecules at deeper levels, infer protection happened
        deep_protected = [
            mol
            for depth, mol in protected_molecules
            if depth > max(d for d, _ in deprotection_reactions)
        ]
        if deep_protected:
            # Infer protection at the deepest level where we see a Boc group
            deepest_protected = max(protected_molecules, key=lambda x: x[0])
            protection_reactions.append(deepest_protected)
            print(f"Inferred Boc protection at depth {deepest_protected[0]}")

    # Still no complete strategy
    if not protection_reactions or not deprotection_reactions:
        print("No complete Boc protection strategy found")
        return False

    # In retrosynthetic analysis:
    # - Early stage = higher depth
    # - Late stage = lower depth
    early_threshold = max_depth * 2 // 3  # Early stage is in the deeper part of the tree
    late_threshold = max_depth // 3  # Late stage is in the shallower part of the tree

    earliest_protection = max(protection_reactions, key=lambda x: x[0])[
        0
    ]  # Highest depth = earliest
    latest_deprotection = min(deprotection_reactions, key=lambda x: x[0])[
        0
    ]  # Lowest depth = latest

    print(
        f"Earliest protection at depth {earliest_protection}, latest deprotection at depth {latest_deprotection}"
    )
    print(f"Early threshold: {early_threshold}, late threshold: {late_threshold}")

    # Check if there are protected intermediates between protection and deprotection
    # In retrosynthetic analysis, deprotection (lower depth) happens before protection (higher depth)
    intermediate_steps = False
    for depth, mol in protected_molecules:
        if latest_deprotection < depth < earliest_protection:
            intermediate_steps = True
            print(f"Found protected intermediate at depth {depth}")
            break

    # A good Boc strategy has:
    # - Early protection (high depth in retrosynthesis)
    # - Late deprotection (low depth in retrosynthesis)
    # - Protected intermediates between them
    is_good_strategy = (
        earliest_protection >= early_threshold
        and latest_deprotection <= late_threshold
        and intermediate_steps
    )

    print(f"Is good Boc protection strategy: {is_good_strategy}")
    return is_good_strategy

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
    Detects a synthesis featuring multiple different nitrogen protecting groups
    (at least 2 different types) being used throughout the synthesis.
    """
    # Track protecting groups and events
    protecting_groups_used = set()
    protection_events = []
    deprotection_events = []

    # List of common nitrogen protecting groups to check
    n_protecting_groups = ["Boc", "Phthalimide", "Tosyl", "Nosyl", "Acetyl", "Trifluoroacetyl"]

    # Import RDKit for substructure matching
    from rdkit import Chem

    # Define Cbz pattern for detection
    cbz_pattern = Chem.MolFromSmiles("C(=O)OCc1ccccc1")

    def has_cbz_group(smiles):
        """Check if a molecule contains a Cbz group"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol and mol.HasSubstructMatch(cbz_pattern):
                return True
        except:
            pass
        return False

    def dfs_traverse(node, depth=0):
        nonlocal protecting_groups_used, protection_events, deprotection_events

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for protecting groups in molecules
            for pg in n_protecting_groups:
                if checker.check_fg(pg, mol_smiles):
                    protecting_groups_used.add(pg)
                    print(f"Found protecting group {pg} in molecule: {mol_smiles}")

            # Check for Cbz group
            if has_cbz_group(mol_smiles):
                protecting_groups_used.add("Cbz")
                print(f"Found Cbz protecting group in molecule: {mol_smiles}")

        elif node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc protection reactions
                if (
                    checker.check_reaction("Boc amine protection", rsmi)
                    or checker.check_reaction("Boc amine protection explicit", rsmi)
                    or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                    or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                    or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                    or checker.check_reaction("Boc amine protection of primary amine", rsmi)
                ):
                    protection_events.append(("Boc", rsmi))
                    protecting_groups_used.add("Boc")
                    print(f"Detected Boc protection reaction: {rsmi}")

                # Check for Boc deprotection reactions
                if (
                    checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                    or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                ):
                    deprotection_events.append(("Boc", rsmi))
                    protecting_groups_used.add("Boc")
                    print(f"Detected Boc deprotection reaction: {rsmi}")

                # Check for Cbz protection/deprotection
                reactants_have_cbz = any(has_cbz_group(r) for r in reactants)
                product_has_cbz = has_cbz_group(product)

                if not reactants_have_cbz and product_has_cbz:
                    protection_events.append(("Cbz", rsmi))
                    protecting_groups_used.add("Cbz")
                    print(f"Detected Cbz protection: {rsmi}")

                if reactants_have_cbz and not product_has_cbz:
                    deprotection_events.append(("Cbz", rsmi))
                    protecting_groups_used.add("Cbz")
                    print(f"Detected Cbz deprotection: {rsmi}")

                # Check for Carboxyl benzyl deprotection reaction
                if checker.check_reaction("Carboxyl benzyl deprotection", rsmi):
                    deprotection_events.append(("Cbz", rsmi))
                    protecting_groups_used.add("Cbz")
                    print(f"Detected Cbz deprotection reaction: {rsmi}")

                # Check for phthalimide protection/deprotection
                reactants_have_phthalimide = any(
                    checker.check_fg("Phthalimide", r) for r in reactants
                )
                product_has_phthalimide = checker.check_fg("Phthalimide", product)

                if not reactants_have_phthalimide and product_has_phthalimide:
                    protection_events.append(("Phthalimide", rsmi))
                    protecting_groups_used.add("Phthalimide")
                    print(f"Detected Phthalimide protection: {rsmi}")

                if reactants_have_phthalimide and not product_has_phthalimide:
                    deprotection_events.append(("Phthalimide", rsmi))
                    protecting_groups_used.add("Phthalimide")
                    print(f"Detected Phthalimide deprotection: {rsmi}")

                if checker.check_reaction("Phthalimide deprotection", rsmi):
                    deprotection_events.append(("Phthalimide", rsmi))
                    protecting_groups_used.add("Phthalimide")
                    print(f"Detected Phthalimide deprotection reaction: {rsmi}")

                # Check for other protection/deprotection events
                for pg in n_protecting_groups:
                    if pg not in ["Boc", "Phthalimide"]:  # Already handled above
                        reactants_have_pg = any(checker.check_fg(pg, r) for r in reactants)
                        product_has_pg = checker.check_fg(pg, product)

                        if not reactants_have_pg and product_has_pg:
                            protection_events.append((pg, rsmi))
                            protecting_groups_used.add(pg)
                            print(f"Detected {pg} protection: {rsmi}")

                        if reactants_have_pg and not product_has_pg:
                            deprotection_events.append((pg, rsmi))
                            protecting_groups_used.add(pg)
                            print(f"Detected {pg} deprotection: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Count total protection/deprotection events
    total_events = len(protection_events) + len(deprotection_events)

    # Check if the strategy is present
    has_multiple_protecting_groups = len(protecting_groups_used) >= 2
    has_multiple_events = total_events >= 2

    print(f"Protecting groups used: {protecting_groups_used}")
    print(f"Protection events: {len(protection_events)}")
    print(f"Deprotection events: {len(deprotection_events)}")
    print(f"Total events: {total_events}")

    # If we have at least 2 different protecting groups in the molecules,
    # consider it a multiple protection strategy even if we don't detect events
    if len(protecting_groups_used) >= 2:
        print(f"Multiple nitrogen protection strategy detected!")
        return True

    # Original condition
    if has_multiple_protecting_groups and has_multiple_events:
        print(f"Multiple nitrogen protection strategy detected!")
        return True

    return False

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
    This function detects if the final step in the synthesis involves
    a deprotection reaction.
    """
    # In retrosynthesis, the root node is the final product
    # The first reaction node is the final step

    # If the route is a molecule node (final product), check its children
    if route["type"] == "mol":
        # Look for the first reaction node (final step)
        for child in route.get("children", []):
            if (
                child["type"] == "reaction"
                and "metadata" in child
                and "rsmi" in child["metadata"]
            ):
                rsmi = child["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking final step reaction: {rsmi}")

                # Check for various deprotection reactions
                if checker.check_reaction("Boc amine deprotection", rsmi):
                    print("Found late-stage deprotection: Boc deprotection")
                    return True

                if checker.check_reaction(
                    "Alcohol deprotection from silyl ethers", rsmi
                ):
                    print("Found late-stage deprotection: Silyl ether deprotection")
                    return True

                if checker.check_reaction(
                    "Ester saponification (methyl deprotection)", rsmi
                ):
                    print("Found late-stage deprotection: Ester saponification")
                    return True

                if checker.check_reaction(
                    "Ester saponification (alkyl deprotection)", rsmi
                ):
                    print("Found late-stage deprotection: Ester saponification")
                    return True

                if checker.check_reaction("COOH ethyl deprotection", rsmi):
                    print("Found late-stage deprotection: Carboxylic acid deprotection")
                    return True

                if checker.check_reaction("Hydroxyl benzyl deprotection", rsmi):
                    print("Found late-stage deprotection: Benzyl deprotection")
                    return True

                if checker.check_reaction("Carboxyl benzyl deprotection", rsmi):
                    print("Found late-stage deprotection: Benzyl deprotection")
                    return True

                if checker.check_reaction(
                    "Cleavage of methoxy ethers to alcohols", rsmi
                ):
                    print("Found late-stage deprotection: Methoxy ether cleavage")
                    return True

                if checker.check_reaction(
                    "Cleavage of alkoxy ethers to alcohols", rsmi
                ):
                    print("Found late-stage deprotection: Alkoxy ether cleavage")
                    return True

                if checker.check_reaction("TMS deprotection from alkyne", rsmi):
                    print("Found late-stage deprotection: TMS alkyne deprotection")
                    return True

                if checker.check_reaction("Tert-butyl deprotection of amine", rsmi):
                    print(
                        "Found late-stage deprotection: Tert-butyl amine deprotection"
                    )
                    return True

                if checker.check_reaction("Phthalimide deprotection", rsmi):
                    print("Found late-stage deprotection: Phthalimide deprotection")
                    return True

                if checker.check_reaction("N-glutarimide deprotection", rsmi):
                    print("Found late-stage deprotection: N-glutarimide deprotection")
                    return True

                if checker.check_reaction("Acetal hydrolysis to aldehyde", rsmi):
                    print("Found late-stage deprotection: Acetal hydrolysis")
                    return True

                if checker.check_reaction("Ketal hydrolysis to ketone", rsmi):
                    print("Found late-stage deprotection: Ketal hydrolysis")
                    return True

                # Check for methoxy to carbonyl conversion by examining atom mapping
                product_mol = Chem.MolFromSmiles(product)
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and product_mol:
                        # Look for specific atom-mapped patterns
                        for atom_idx in range(
                            1, 100
                        ):  # Check reasonable range of atom indices
                            methoxy_pattern = f"[O:{atom_idx}]"
                            carbonyl_pattern = f"=[O:{atom_idx}]"
                            if (
                                methoxy_pattern in reactant
                                and carbonyl_pattern in product
                            ):
                                print(
                                    f"Found late-stage deprotection: Methoxy to carbonyl conversion (atom-mapped at position {atom_idx})"
                                )
                                return True

                # Check for methoxy to carbonyl conversion (common deprotection)
                for reactant in reactants:
                    # Check for methoxy group in reactant
                    if checker.check_fg("Ether", reactant):
                        # Check for carbonyl in product that wasn't in reactant
                        if (
                            checker.check_fg("Ketone", product)
                            and not checker.check_fg("Ketone", reactant)
                        ) or (
                            checker.check_fg("Aldehyde", product)
                            and not checker.check_fg("Aldehyde", reactant)
                        ):
                            print(
                                "Found late-stage deprotection: Methoxy to carbonyl conversion"
                            )
                            return True

                # Check for ester to amide conversion (not strictly a deprotection)
                for reactant in reactants:
                    if (
                        checker.check_fg("Ester", reactant)
                        and not checker.check_fg("Ester", product)
                        and checker.check_fg("Primary amide", product)
                    ):
                        print(
                            "Found late-stage deprotection: Ester to amide conversion"
                        )
                        return True

                # Check for protective group removal
                for reactant in reactants:
                    # Check for silyl protective group removal
                    if checker.check_fg(
                        "Silyl protective group", reactant
                    ) and not checker.check_fg("Silyl protective group", product):
                        print(
                            "Found late-stage deprotection: Silyl protective group removal"
                        )
                        return True

                    # Check for TMS ether protective group removal
                    if checker.check_fg(
                        "TMS ether protective group", reactant
                    ) and not checker.check_fg("TMS ether protective group", product):
                        print(
                            "Found late-stage deprotection: TMS ether protective group removal"
                        )
                        return True

                    # Check for Boc removal
                    if checker.check_fg("Boc", reactant) and not checker.check_fg(
                        "Boc", product
                    ):
                        print("Found late-stage deprotection: Boc removal")
                        return True

                    # Check for acetal/ketal deprotection
                    if checker.check_fg(
                        "Acetal/Ketal", reactant
                    ) and not checker.check_fg("Acetal/Ketal", product):
                        print(
                            "Found late-stage deprotection: Acetal/Ketal deprotection"
                        )
                        return True

                # Check for general deprotection patterns
                if any(
                    term in rsmi.lower()
                    for term in ["deprotection", "deprotect", "hydrolysis", "cleavage"]
                ):
                    print(
                        "Found late-stage deprotection: Deprotection mentioned in reaction SMILES"
                    )
                    return True

    return False

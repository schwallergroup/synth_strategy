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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects if the synthesis involves late-stage aromatic bromination.
    """
    late_bromination_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal late_bromination_detected

        if node["type"] == "reaction" and depth <= 3 and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Check for aromatic bromination reaction using the provided checker
            if checker.check_reaction("Aromatic bromination", rsmi):
                print(f"Late-stage aromatic bromination detected at depth {depth}")
                late_bromination_detected = True
                return

            # If the direct check failed, analyze the reaction more carefully
            try:
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Split reactants if there are multiple
                reactants = reactants_part.split(".")

                # Convert to RDKit molecules
                product_mol = Chem.MolFromSmiles(product_part)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]

                if product_mol and reactant_mols:
                    # Check if product contains aromatic-Br bonds
                    ar_br_pattern = Chem.MolFromSmarts("c-[Br]")
                    product_ar_br_count = len(product_mol.GetSubstructMatches(ar_br_pattern))

                    if product_ar_br_count > 0:
                        print(f"Product contains {product_ar_br_count} aromatic-Br bonds")

                        # Check if any reactant contains aromatic rings
                        for i, reactant_mol in enumerate(reactant_mols):
                            if reactant_mol:
                                reactant_ar_br_count = len(
                                    reactant_mol.GetSubstructMatches(ar_br_pattern)
                                )
                                print(
                                    f"Reactant {i} contains {reactant_ar_br_count} aromatic-Br bonds"
                                )

                                # If product has more aromatic-Br bonds than this reactant,
                                # and the reactant has aromatic rings, it might be bromination
                                if product_ar_br_count > reactant_ar_br_count:
                                    # Check if reactant has aromatic rings
                                    aromatic_pattern = Chem.MolFromSmarts("c")
                                    if reactant_mol.HasSubstructMatch(aromatic_pattern):
                                        # Check if any reactant contains bromine source
                                        br_sources = ["Br", "CBr", "HBr", "BBr", "NBr"]
                                        br_source_present = any(
                                            any(br_src in r for br_src in br_sources)
                                            for r in reactants
                                        )

                                        if br_source_present:
                                            print(
                                                f"Late-stage aromatic bromination detected through detailed analysis at depth {depth}"
                                            )
                                            late_bromination_detected = True
                                            return

                # Additional check for specific bromination reactions
                if (
                    any("c" in r and "Br" in r for r in reactants)
                    and "c" in product_part
                    and "Br" in product_part
                ):
                    # Check if this is a bromination reaction by name
                    if (
                        checker.check_reaction("Wohl-Ziegler bromination benzyl primary", rsmi)
                        or checker.check_reaction("Wohl-Ziegler bromination benzyl secondary", rsmi)
                        or checker.check_reaction("Wohl-Ziegler bromination benzyl tertiary", rsmi)
                        or checker.check_reaction("Wohl-Ziegler bromination allyl primary", rsmi)
                        or checker.check_reaction("Wohl-Ziegler bromination allyl secondary", rsmi)
                        or checker.check_reaction("Wohl-Ziegler bromination allyl tertiary", rsmi)
                    ):
                        print(f"Late-stage Wohl-Ziegler bromination detected at depth {depth}")
                        late_bromination_detected = True
                        return
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {late_bromination_detected}")

    return late_bromination_detected

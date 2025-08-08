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
    This function detects a "build-up" synthetic strategy where functional groups
    are sequentially added to a core structure that remains intact throughout.
    """
    # Track the number of functional group additions
    functional_group_additions = 0
    core_modified = False

    # List of functional groups to check for additions
    functional_groups = [
        "Aromatic halide",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Nitro group",
        "Ester",
        "Nitrile",
        "Carboxylic acid",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Alcohol",
        "Ether",
        "Ketone",
        "Aldehyde",
        "Amide",
        "Sulfonamide",
        "Azide",
        "Alkyne",
        "Boronic acid",
        "Boronic ester",
    ]

    # List of reaction types that typically add functional groups
    fg_addition_reactions = [
        "Aromatic bromination",
        "Aromatic chlorination",
        "Aromatic fluorination",
        "Aromatic iodination",
        "Aromatic nitration",
        "Esterification",
        "Acylation",
        "Alkylation",
        "Suzuki coupling",
        "Sonogashira",
    ]

    # Core structures to check
    core_structures = [
        "benzene",
        "naphthalene",
        "anthracene",
        "pyridine",
        "pyrimidine",
        "furan",
        "thiophene",
        "pyrrole",
        "indole",
        "quinoline",
        "isoquinoline",
    ]

    # Track which core structures are present in the final product
    final_product_cores = set()

    def identify_cores(smiles):
        """Identify core structures in a molecule"""
        cores_found = set()
        for core in core_structures:
            if checker.check_ring(core, smiles):
                cores_found.add(core)
        return cores_found

    def dfs_traverse(node, depth=0):
        nonlocal functional_group_additions, core_modified, final_product_cores

        # If this is the final product (depth 0), identify its core structures
        if depth == 0 and node["type"] == "mol":
            final_product_cores = identify_cores(node["smiles"])
            print(f"Final product cores: {final_product_cores}")

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for functional group additions by comparing product and reactants
            product_fgs = set()
            reactant_fgs = set()

            for fg in functional_groups:
                if checker.check_fg(fg, product):
                    product_fgs.add(fg)

                for reactant in reactants:
                    if checker.check_fg(fg, reactant):
                        reactant_fgs.add(fg)

            # Check if any new functional groups were added
            new_fgs = product_fgs - reactant_fgs
            if new_fgs:
                functional_group_additions += 1
                print(f"Found functional group addition: {new_fgs}")

            # Also check if the reaction type is known to add functional groups
            for rxn_type in fg_addition_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    print(f"Detected functional group addition reaction: {rxn_type}")
                    functional_group_additions += 1
                    break

            # Check if core structures are preserved
            product_cores = identify_cores(product)
            reactant_cores = set()
            for reactant in reactants:
                reactant_cores.update(identify_cores(reactant))

            # If we're in an early stage (high depth) and core structures change significantly,
            # this might be core formation rather than modification
            if depth > 2:
                # Early stage - core formation is expected
                pass
            else:
                # Late stage - core should be preserved
                if final_product_cores and not any(
                    core in product_cores for core in final_product_cores
                ):
                    core_modified = True
                    print(f"Core structure was modified at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if we have multiple functional group additions without core modification
    return functional_group_additions >= 2 and not core_modified

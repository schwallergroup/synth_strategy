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
    This function detects a functional group interconversion sequence:
    ketone reduction followed by alcohol protection.
    """
    # Track reactions in sequence
    reduction_steps = []
    protection_steps = []

    def dfs_traverse(node, depth=0, path=[]):
        nonlocal reduction_steps, protection_steps

        current_path = path + [node]

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    # Check for various reduction reactions
                    is_reduction = False

                    # Check for ketone reduction
                    if checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi):
                        is_reduction = True
                    # Check for aldehyde reduction
                    elif checker.check_reaction(
                        "Reduction of aldehydes and ketones to alcohols", rsmi
                    ):
                        is_reduction = True
                    # Check for ester reduction
                    elif checker.check_reaction("Reduction of ester to primary alcohol", rsmi):
                        is_reduction = True
                    # Check for carboxylic acid reduction
                    elif checker.check_reaction(
                        "Reduction of carboxylic acid to primary alcohol", rsmi
                    ):
                        is_reduction = True

                    if is_reduction:
                        # Verify functional groups
                        reactant_has_carbonyl = any(
                            checker.check_fg("Ketone", reactant)
                            or checker.check_fg("Aldehyde", reactant)
                            or checker.check_fg("Ester", reactant)
                            or checker.check_fg("Carboxylic acid", reactant)
                            for reactant in reactants
                        )

                        product_has_alcohol = (
                            checker.check_fg("Secondary alcohol", product)
                            or checker.check_fg("Primary alcohol", product)
                            or checker.check_fg("Tertiary alcohol", product)
                        )

                        if reactant_has_carbonyl and product_has_alcohol:
                            reduction_steps.append((depth, node, product))
                            print(
                                f"Carbonyl reduction detected at depth {depth}, product: {product}"
                            )

                    # Check for alcohol protection
                    alcohol_reactant = None
                    for reactant in reactants:
                        if (
                            checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                        ):
                            alcohol_reactant = reactant
                            break

                    # Check for various protection reactions
                    is_protection = False
                    if alcohol_reactant:
                        # Check for silyl ether protection
                        if checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                            is_protection = True
                        # Check for acetalization
                        elif checker.check_reaction("Aldehyde or ketone acetalization", rsmi):
                            is_protection = True
                        # Check for esterification
                        elif checker.check_reaction("Esterification of Carboxylic Acids", rsmi):
                            is_protection = True
                        # Check for TMS protection
                        elif checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                            is_protection = True
                        # Check for ether formation
                        elif checker.check_reaction("Williamson Ether Synthesis", rsmi):
                            is_protection = True
                        # Check for acetal formation
                        elif checker.check_reaction("Diol acetalization", rsmi):
                            is_protection = True
                        # Check for other protection reactions where alcohol is converted
                        elif (
                            not checker.check_fg("Secondary alcohol", product)
                            and not checker.check_fg("Primary alcohol", product)
                            and not checker.check_fg("Tertiary alcohol", product)
                        ):
                            # The alcohol is gone, likely protected
                            is_protection = True

                    if is_protection and alcohol_reactant:
                        protection_steps.append((depth, node, alcohol_reactant))
                        print(
                            f"Alcohol protection detected at depth {depth}, reactant: {alcohol_reactant}"
                        )
                except Exception as e:
                    print(f"Error in processing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

    # Start traversal from root
    dfs_traverse(route)

    print(
        f"Found {len(reduction_steps)} reduction steps and {len(protection_steps)} protection steps"
    )

    # Check if reduction is followed by protection in the synthetic sequence
    for red_depth, red_node, red_product in reduction_steps:
        for prot_depth, prot_node, prot_reactant in protection_steps:
            # Protection should occur after reduction (higher depth in retrosynthesis)
            if prot_depth > red_depth:
                print(
                    f"Checking reduction at depth {red_depth} and protection at depth {prot_depth}"
                )

                # Check if there's a path from reduction to protection
                def is_in_path(start_node, target_node, visited=None):
                    if visited is None:
                        visited = set()

                    if start_node == target_node:
                        return True

                    visited.add(start_node)

                    for child in start_node.get("children", []):
                        if child not in visited and is_in_path(child, target_node, visited):
                            return True

                    return False

                # Check if protection node is in the path from reduction node
                if is_in_path(red_node, prot_node):
                    print("Protection node is in the path from reduction node")

                    # Compare the molecules
                    red_mol = Chem.MolFromSmiles(red_product)
                    prot_mol = Chem.MolFromSmiles(prot_reactant)

                    if red_mol and prot_mol:
                        # Use MCS to compare molecules
                        mcs = rdFMCS.FindMCS(
                            [red_mol, prot_mol],
                            bondCompare=rdFMCS.BondCompare.CompareOrder,
                            atomCompare=rdFMCS.AtomCompare.CompareElements,
                            ringMatchesRingOnly=True,
                            completeRingsOnly=True,
                        )

                        red_atom_count = red_mol.GetNumAtoms()
                        prot_atom_count = prot_mol.GetNumAtoms()
                        mcs_atom_count = mcs.numAtoms

                        # Calculate similarity as ratio of MCS atoms to smaller molecule's atoms
                        similarity = mcs_atom_count / min(red_atom_count, prot_atom_count)
                        print(
                            f"Molecule similarity: {similarity:.2f} (MCS atoms: {mcs_atom_count}, Red atoms: {red_atom_count}, Prot atoms: {prot_atom_count})"
                        )

                        # Check if the alcohol group from reduction is the same one being protected
                        # If molecules are similar enough (70% or more atoms in common)
                        if similarity >= 0.7:
                            # Check if the alcohol functional group is present in both molecules
                            red_has_alcohol = (
                                checker.check_fg("Secondary alcohol", red_product)
                                or checker.check_fg("Primary alcohol", red_product)
                                or checker.check_fg("Tertiary alcohol", red_product)
                            )

                            prot_has_alcohol = (
                                checker.check_fg("Secondary alcohol", prot_reactant)
                                or checker.check_fg("Primary alcohol", prot_reactant)
                                or checker.check_fg("Tertiary alcohol", prot_reactant)
                            )

                            if red_has_alcohol and prot_has_alcohol:
                                print(
                                    f"Reduction-protection sequence detected: reduction at depth {red_depth}, protection at depth {prot_depth}"
                                )
                                return True

    # If we have protection steps but no reduction steps, check if any protection
    # involves a molecule that could have come from a reduction
    if len(reduction_steps) == 0 and len(protection_steps) > 0:
        for prot_depth, prot_node, prot_reactant in protection_steps:
            # Check if the protected alcohol could have come from a reduction
            if (
                checker.check_fg("Secondary alcohol", prot_reactant)
                or checker.check_fg("Primary alcohol", prot_reactant)
                or checker.check_fg("Tertiary alcohol", prot_reactant)
            ):

                # Look for parent nodes that might be reductions but weren't caught
                def check_parent_for_reduction(node, target, depth=0, path=[]):
                    if node == target:
                        return False

                    if node["type"] == "reaction":
                        if "rsmi" in node.get("metadata", {}):
                            rsmi = node["metadata"]["rsmi"]
                            product = rsmi.split(">")[-1]

                            # Check if product has an alcohol group
                            has_alcohol = (
                                checker.check_fg("Secondary alcohol", product)
                                or checker.check_fg("Primary alcohol", product)
                                or checker.check_fg("Tertiary alcohol", product)
                            )

                            if has_alcohol:
                                # This might be a reduction that wasn't caught
                                print(
                                    f"Potential uncaught reduction at depth {depth}, product: {product}"
                                )

                                # Compare with protection reactant
                                prod_mol = Chem.MolFromSmiles(product)
                                prot_mol = Chem.MolFromSmiles(prot_reactant)

                                if prod_mol and prot_mol:
                                    mcs = rdFMCS.FindMCS(
                                        [prod_mol, prot_mol],
                                        bondCompare=rdFMCS.BondCompare.CompareOrder,
                                        atomCompare=rdFMCS.AtomCompare.CompareElements,
                                    )

                                    similarity = mcs.numAtoms / min(
                                        prod_mol.GetNumAtoms(), prot_mol.GetNumAtoms()
                                    )
                                    print(f"Similarity with protection reactant: {similarity:.2f}")

                                    if similarity >= 0.7:
                                        print(f"Potential reduction-protection sequence detected")
                                        return True

                    for child in node.get("children", []):
                        if check_parent_for_reduction(child, target, depth + 1, path + [node]):
                            return True

                    return False

                if check_parent_for_reduction(route, prot_node):
                    return True

    return False

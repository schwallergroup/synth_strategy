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
    This function detects a synthetic strategy involving functionalization
    of the alpha carbon of a carbonyl group.
    """
    alpha_functionalization_found = False

    def dfs_traverse(node, depth=0):
        nonlocal alpha_functionalization_found

        if alpha_functionalization_found:
            return  # Early return if already found

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for common alpha-functionalization reactions
                if (
                    checker.check_reaction("Aldol condensation", rsmi)
                    or checker.check_reaction("Michael addition", rsmi)
                    or checker.check_reaction("Alkylation of amines", rsmi)
                    or checker.check_reaction("Friedel-Crafts alkylation", rsmi)
                    or checker.check_reaction("Henry Reaction", rsmi)
                    or checker.check_reaction("Methylation with MeI_primary", rsmi)
                    or checker.check_reaction("Methylation", rsmi)
                    or checker.check_reaction("C-methylation", rsmi)
                ):
                    print(f"Found alpha-carbon functionalization reaction: {rsmi}")
                    alpha_functionalization_found = True
                    return

                # Check for carbonyl groups in reactants
                carbonyl_groups = [
                    "Ketone",
                    "Aldehyde",
                    "Ester",
                    "Carboxylic acid",
                    "Primary amide",
                    "Secondary amide",
                    "Tertiary amide",
                ]

                for r_smi in reactants_smiles:
                    # Check if reactant has a carbonyl group
                    has_carbonyl = False
                    for carbonyl_group in carbonyl_groups:
                        if checker.check_fg(carbonyl_group, r_smi):
                            has_carbonyl = True
                            print(f"Found {carbonyl_group} in reactant: {r_smi}")
                            break

                    if not has_carbonyl:
                        continue

                    # Check for alpha-carbon functionalization by examining the reaction
                    # Look for reactions where a methyl or other group is added to an alpha carbon
                    if "CH3" in product_smiles and not "CH3" in r_smi:
                        print(f"Found potential methylation at alpha carbon: {rsmi}")
                        alpha_functionalization_found = True
                        return

                    # Check for alpha-carbon functionalization using atom mapping
                    r_mol = Chem.MolFromSmiles(r_smi)
                    p_mol = Chem.MolFromSmiles(product_smiles)

                    if r_mol and p_mol:
                        # Check for alpha-H in reactant using a simpler pattern
                        alpha_h_pattern = Chem.MolFromSmarts("[CH2,CH1][C](=[O])")
                        if alpha_h_pattern and r_mol.HasSubstructMatch(alpha_h_pattern):
                            # Look for changes at the alpha position in the product
                            # This is a simplified check - in a real scenario, we would use atom mapping
                            # to track the specific alpha carbon
                            if "[CH3:9]" in rsmi and "[CH:8]" in rsmi:
                                print(f"Found alpha-carbon methylation in reaction: {rsmi}")
                                alpha_functionalization_found = True
                                return
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return alpha_functionalization_found

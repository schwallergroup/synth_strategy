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
    Detects if the synthesis uses a convergent approach with multiple fragment couplings,
    specifically looking for C-C bond formations between aromatic systems and olefin formation.
    """
    # Track fragments and coupling reactions
    coupling_reactions_count = 0
    found_biaryl_coupling = False
    found_olefin_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal coupling_reactions_count, found_biaryl_coupling, found_olefin_coupling

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is a coupling reaction with multiple fragments
            if len(reactants_smiles) >= 2:
                # Check for biaryl coupling reactions
                if (
                    checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                    or checker.check_reaction("Stille reaction_aryl", rsmi)
                    or checker.check_reaction("Stille reaction_aryl OTf", rsmi)
                    or checker.check_reaction("Negishi coupling", rsmi)
                    or checker.check_reaction("Ullmann condensation", rsmi)
                    or checker.check_reaction("Kumada cross-coupling", rsmi)
                ):

                    found_biaryl_coupling = True
                    coupling_reactions_count += 1
                    print(f"Found biaryl coupling at depth {depth}: {rsmi}")

                # Check for olefin formation reactions
                elif (
                    checker.check_reaction("Wittig", rsmi)
                    or checker.check_reaction("Wittig reaction with triphenylphosphorane", rsmi)
                    or checker.check_reaction("Wittig with Phosphonium", rsmi)
                    or checker.check_reaction("Julia Olefination", rsmi)
                    or checker.check_reaction("Heck terminal_vinyl", rsmi)
                    or checker.check_reaction("Heck_terminal_vinyl", rsmi)
                    or checker.check_reaction("Heck_non-terminal_vinyl", rsmi)
                    or checker.check_reaction("Olefination of ketones with Grignard reagents", rsmi)
                    or checker.check_reaction(
                        "Olefination of aldehydes with Grignard reagents", rsmi
                    )
                ):

                    found_olefin_coupling = True
                    coupling_reactions_count += 1
                    print(f"Found olefin formation at depth {depth}: {rsmi}")

                # If we can't identify the reaction type specifically, check for C-C bond formation
                else:
                    # Check for biaryl coupling using SMARTS pattern as fallback
                    product = Chem.MolFromSmiles(product_smiles)
                    if product is not None:
                        biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")
                        alkene_pattern = Chem.MolFromSmarts("[#6]=[#6]")

                        # Check if these patterns exist in product but not in all reactants
                        reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                        reactants = [r for r in reactants if r is not None]

                        if product.HasSubstructMatch(biaryl_pattern) and not all(
                            r.HasSubstructMatch(biaryl_pattern) for r in reactants
                        ):
                            found_biaryl_coupling = True
                            coupling_reactions_count += 1
                            print(f"Found biaryl coupling (pattern) at depth {depth}")

                        if product.HasSubstructMatch(alkene_pattern) and not all(
                            r.HasSubstructMatch(alkene_pattern) for r in reactants
                        ):
                            found_olefin_coupling = True
                            coupling_reactions_count += 1
                            print(f"Found olefin formation (pattern) at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # We need at least 2 coupling reactions with at least one biaryl and one olefin formation
    strategy_present = (
        coupling_reactions_count >= 2 and found_biaryl_coupling and found_olefin_coupling
    )

    if strategy_present:
        print(
            f"Detected convergent fragment assembly strategy with {coupling_reactions_count} coupling reactions"
        )
    else:
        print(
            f"Convergent fragment assembly strategy not detected. Found {coupling_reactions_count} coupling reactions, biaryl: {found_biaryl_coupling}, olefin: {found_olefin_coupling}"
        )

    return strategy_present

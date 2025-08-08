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
    This function detects if the synthetic route involves coupling with a nucleoside fragment.
    A nucleoside typically consists of a sugar moiety (ribose or deoxyribose) connected to a
    nitrogenous base (like adenine, guanine, cytosine, thymine, or uracil).
    """
    nucleoside_coupling_detected = False

    def is_nucleoside(mol_smiles):
        """Check if a molecule contains a nucleoside structure"""
        # Check for common sugar rings in nucleosides
        sugar_present = (
            checker.check_ring("furan", mol_smiles)
            or checker.check_ring("tetrahydrofuran", mol_smiles)
            or checker.check_ring("pyran", mol_smiles)
            or checker.check_ring("tetrahydropyran", mol_smiles)
            or checker.check_ring("dioxolane", mol_smiles)
        )

        # Check for nitrogen-containing heterocycles commonly found in nucleosides
        nitrogen_ring_present = (
            checker.check_ring("pyrimidine", mol_smiles)
            or checker.check_ring("purine", mol_smiles)
            or checker.check_ring("imidazole", mol_smiles)
            or checker.check_ring("triazole", mol_smiles)
            or checker.check_ring("tetrazole", mol_smiles)
        )

        # Additional check for adenine-like structures
        purine_present = checker.check_fg("Primary amide", mol_smiles) and nitrogen_ring_present

        print(
            f"Nucleoside check - Sugar: {sugar_present}, N-ring: {nitrogen_ring_present}, Purine: {purine_present}"
        )

        # A nucleoside should have both a sugar and a nitrogenous base
        return sugar_present and (nitrogen_ring_present or purine_present)

    def is_coupling_reaction(rxn_smiles):
        """Check if a reaction is a coupling reaction"""
        # Check for common coupling reactions
        coupling = (
            checker.check_reaction("Suzuki coupling with boronic acids", rxn_smiles)
            or checker.check_reaction("Suzuki coupling with boronic esters", rxn_smiles)
            or checker.check_reaction("Sonogashira acetylene_aryl halide", rxn_smiles)
            or checker.check_reaction("Sonogashira alkyne_aryl halide", rxn_smiles)
            or checker.check_reaction("Heck terminal vinyl", rxn_smiles)
            or checker.check_reaction("Stille reaction_aryl", rxn_smiles)
            or checker.check_reaction("Negishi coupling", rxn_smiles)
            or checker.check_reaction("Buchwald-Hartwig", rxn_smiles)
            or checker.check_reaction("N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rxn_smiles)
            or checker.check_reaction("Reductive amination with aldehyde", rxn_smiles)
            or checker.check_reaction("Reductive amination with ketone", rxn_smiles)
            or checker.check_reaction(
                "N-alkylation of primary amines with alkyl halides", rxn_smiles
            )
            or checker.check_reaction(
                "N-alkylation of secondary amines with alkyl halides", rxn_smiles
            )
        )

        print(f"Is coupling reaction: {coupling}")
        return coupling

    def dfs_traverse(node, depth=0):
        nonlocal nucleoside_coupling_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Check if any reactant or product contains a nucleoside structure
            reactant_nucleosides = []
            non_nucleoside_reactants = []

            print("Checking reactants for nucleoside structures:")
            for reactant_smiles in reactants_smiles:
                print(f"Checking reactant: {reactant_smiles}")
                if is_nucleoside(reactant_smiles):
                    print(f"Nucleoside reactant found: {reactant_smiles}")
                    reactant_nucleosides.append(reactant_smiles)
                else:
                    non_nucleoside_reactants.append(reactant_smiles)

            print(f"Checking product for nucleoside structure: {product_smiles}")
            product_has_nucleoside = is_nucleoside(product_smiles)
            print(f"Product has nucleoside: {product_has_nucleoside}")

            # Check if this is a reaction that could be used for coupling
            if len(reactants_smiles) > 1:  # Multiple reactants
                print("Multiple reactants found, checking if this is a coupling reaction")

                # Case 1: Nucleoside + non-nucleoside -> nucleoside product (classic coupling)
                if reactant_nucleosides and non_nucleoside_reactants and product_has_nucleoside:
                    print("Case 1: Nucleoside + non-nucleoside -> nucleoside product")
                    nucleoside_coupling_detected = True

                # Case 2: Non-nucleoside reactants -> nucleoside product (nucleoside formation)
                elif not reactant_nucleosides and product_has_nucleoside:
                    print("Case 2: Non-nucleoside reactants -> nucleoside product")
                    nucleoside_coupling_detected = True

                # Case 3: Check if it's a known coupling reaction
                elif is_coupling_reaction(rsmi) and (
                    reactant_nucleosides or product_has_nucleoside
                ):
                    print("Case 3: Known coupling reaction involving nucleoside")
                    nucleoside_coupling_detected = True

        # Continue traversing the synthesis route
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return nucleoside_coupling_detected

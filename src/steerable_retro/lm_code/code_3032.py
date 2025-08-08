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
    This function detects formylation of a heterocycle.
    """
    found_formylation = False

    # List of common heterocycles to check
    heterocycles = [
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "furan",
        "thiophene",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
    ]

    # Formylating agents and related reagents
    formylating_agents = [
        "CN(C)C=O",  # DMF
        "OC=O",  # Formic acid
        "ClC=O",  # Formyl chloride
        "O=CH",  # Formaldehyde
        "HC(=O)O",  # Formate
        "C(=O)H",  # Aldehyde pattern
        "O=CNO",  # Formamide N-oxide
        "HCOOH",  # Formic acid
        "HCOONa",  # Sodium formate
        "CO",  # Carbon monoxide (for metal-catalyzed formylation)
    ]

    def dfs_traverse(node):
        nonlocal found_formylation

        if found_formylation:
            return

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Check if product contains an aldehyde group
                if checker.check_fg("Aldehyde", product):
                    print(f"Found aldehyde in product: {product}")

                    # Check if product contains a heterocycle
                    product_heterocycle = None
                    for ring in heterocycles:
                        if checker.check_ring(ring, product):
                            product_heterocycle = ring
                            print(f"Found heterocycle {ring} in product")
                            break

                    if product_heterocycle:
                        # Check if aldehyde was added in this reaction (not present in reactants)
                        aldehyde_in_reactants = False
                        for reactant in reactants:
                            if checker.check_fg("Aldehyde", reactant):
                                # If a reactant already has an aldehyde, it might not be formylation
                                aldehyde_in_reactants = True
                                print(f"Found aldehyde in reactant: {reactant}")
                                break

                        # Check for formylating agents in reactants
                        formylating_agent_present = False
                        for reactant in reactants:
                            # Check for known formylating agents
                            if any(agent in reactant for agent in formylating_agents):
                                formylating_agent_present = True
                                print(f"Found formylating agent in reactant: {reactant}")
                                break
                            # Check for POCl3 (often used with DMF in Vilsmeier-Haack reaction)
                            elif (
                                "POCl3" in reactant
                                or "P(=O)(Cl)Cl" in reactant
                                or "PCl3" in reactant
                            ):
                                formylating_agent_present = True
                                print(f"Found POCl3 or similar reagent in reactant: {reactant}")
                                break

                        # Check for specific formylation reactions
                        is_formylation_reaction = False

                        # Check if this is a known formylation reaction type
                        if checker.check_reaction("Aromatic hydroxylation", rsmi):
                            print(
                                "Detected aromatic hydroxylation reaction (can be related to formylation)"
                            )
                            is_formylation_reaction = True

                        # If we have a heterocycle with a new aldehyde and a formylating agent, it's likely formylation
                        if not aldehyde_in_reactants and (
                            formylating_agent_present or is_formylation_reaction
                        ):
                            print("Found formylation of heterocycle")
                            found_formylation = True
                            return

                        # Special case: Check for Vilsmeier-Haack reaction pattern
                        # (DMF + POCl3 + heterocycle → formylated heterocycle)
                        dmf_present = any("CN(C)C=O" in reactant for reactant in reactants)
                        pocl3_present = any(
                            ("POCl3" in reactant or "P(=O)(Cl)Cl" in reactant)
                            for reactant in reactants
                        )

                        if dmf_present and pocl3_present and not aldehyde_in_reactants:
                            print("Detected Vilsmeier-Haack formylation pattern")
                            found_formylation = True
                            return

                        # If we've reached this point, check if the reaction is simply adding an aldehyde
                        # to a heterocycle, which is a common formylation pattern
                        if not aldehyde_in_reactants:
                            print("Detected addition of aldehyde to heterocycle")
                            found_formylation = True
                            return
            except Exception as e:
                print(f"Error in processing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_formylation

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
    Detects if the synthesis route includes a late-stage formylation
    (addition of CHO group in the final or penultimate step)
    """
    formylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal formylation_detected

        if node["type"] == "reaction" and depth <= 2:  # Check final, penultimate, and one more step
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if any reactant is a formylating agent
                formyl_agent = False
                for reactant in reactants:
                    # Common formylating agents
                    if checker.check_fg("Formaldehyde", reactant) or (
                        ("C=O" in reactant or "CH=O" in reactant) and len(reactant) < 15
                    ):  # Simple formyl sources
                        formyl_agent = True
                        print(f"Potential formylating agent detected: {reactant}")
                        break

                # Check if product contains aldehyde group but reactants don't all have it
                product_has_aldehyde = checker.check_fg("Aldehyde", product)
                reactants_with_aldehyde = [r for r in reactants if checker.check_fg("Aldehyde", r)]

                # Check for specific formylation reactions
                formylation_reaction = False

                # Common formylation reactions
                if checker.check_reaction("Friedel-Crafts acylation", rsmi):
                    # For Friedel-Crafts, verify it's specifically formylation (adding CHO)
                    if formyl_agent or any("C(=O)H" in r for r in reactants):
                        print(f"Friedel-Crafts formylation detected at depth {depth}")
                        formylation_reaction = True

                # Check for Vilsmeier-Haack formylation (DMF reagent)
                if any(r.find("CN(C)C=O") >= 0 or r.find("CN(C)CH=O") >= 0 for r in reactants):
                    print(f"Vilsmeier-Haack formylation detected at depth {depth}")
                    formylation_reaction = True

                # Check for direct formylation reactions
                if "formylation" in rsmi.lower() or "formyl" in rsmi.lower():
                    print(f"Explicit formylation reaction at depth {depth}")
                    formylation_reaction = True

                # Detect formylation by checking if product has new aldehyde group
                new_aldehyde_formed = False
                if product_has_aldehyde and len(reactants_with_aldehyde) < len(reactants):
                    # At least one reactant doesn't have aldehyde, suggesting new formation
                    print(f"New aldehyde group detected in product at depth {depth}")
                    new_aldehyde_formed = True

                # Determine if this is a formylation reaction
                if formylation_reaction or new_aldehyde_formed:
                    print(f"Late-stage formylation detected at depth {depth}")
                    formylation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return formylation_detected

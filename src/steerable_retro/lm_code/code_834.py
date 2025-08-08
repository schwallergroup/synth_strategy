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
    Detects if the synthesis uses a late-stage Suzuki coupling (depth 0 or 1)
    for final diversification.
    """
    suzuki_at_late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_at_late_stage

        if node["type"] == "reaction" and depth <= 1:
            # Get reaction SMILES
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                print(f"No reaction SMILES found at depth {depth}")
                return

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Check for Suzuki coupling reactions using the checker function
            suzuki_reaction_types = [
                "Suzuki coupling with boronic acids",
                "Suzuki coupling with boronic acids OTf",
                "Suzuki coupling with sulfonic esters",
                "Suzuki coupling with boronic esters OTf",
                "Suzuki coupling with boronic esters",
                "Suzuki",
            ]

            # First try with the checker function
            for reaction_type in suzuki_reaction_types:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Found Suzuki coupling ({reaction_type}) at depth {depth}")
                    suzuki_at_late_stage = True
                    return

            # If checker fails, try to manually identify Suzuki coupling characteristics
            try:
                # Split reaction SMILES to get reactants and product
                reactants_str = rsmi.split(">")[0]
                reagents_str = rsmi.split(">")[1] if len(rsmi.split(">")) > 2 else ""
                product_str = rsmi.split(">")[-1]

                # Check for key components of Suzuki coupling
                reactants = reactants_str.split(".")
                reagents = reagents_str.split(".") if reagents_str else []

                # Check for boronic acid/ester in reactants
                has_boronic = False
                has_aryl_halide = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Check for boronic acid/ester
                        if "B" in reactant and (
                            "OB" in reactant or "B(" in reactant or "BO" in reactant
                        ):
                            print(f"Found boronic acid/ester: {reactant}")
                            has_boronic = True

                        # Check for aryl halide
                        if any(x in reactant for x in ["Br", "I", "Cl"]) and any(
                            x in reactant for x in ["c", "C:"] + list("cdefghijklmno")
                        ):
                            print(f"Found aryl halide: {reactant}")
                            has_aryl_halide = True

                # Check for palladium catalyst in reagents
                has_pd = any("Pd" in reagent for reagent in reagents)

                # If we have both boronic compound and aryl halide, and possibly Pd, it's likely a Suzuki coupling
                if has_boronic and has_aryl_halide:
                    print(
                        f"Manual detection: Found Suzuki coupling characteristics at depth {depth}"
                    )
                    suzuki_at_late_stage = True
                    return

            except Exception as e:
                print(f"Error in manual Suzuki detection: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Late stage Suzuki coupling detected: {suzuki_at_late_stage}")
    return suzuki_at_late_stage

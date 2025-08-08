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
    This function detects a synthetic strategy involving late-stage arylation
    with a heterocyclic fragment (pyridine, furan, thiophene, etc.).
    """
    found_late_arylation = False

    # List of heterocycles to check
    heterocycles = [
        "pyridine",
        "furan",
        "thiophene",
        "pyrrole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrazole",
        "isoxazole",
        "isothiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
    ]

    # List of arylation reaction types
    arylation_reactions = [
        "Suzuki",
        "Negishi",
        "Buchwald-Hartwig",
        "N-arylation",
        "Heck",
        "Stille",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
        "Ullmann-Goldberg Substitution amine",
        "Ullmann-Goldberg Substitution thiol",
        "Ullmann-Goldberg Substitution aryl alcohol",
        "Ullmann condensation",
        "Chan-Lam amine",
        "Chan-Lam alcohol",
        "Chan-Lam etherification",
    ]

    # Define aryl functional groups
    aryl_groups = ["Phenol", "Aniline", "Aromatic halide"]

    def has_aryl_group(smiles):
        """Check if a molecule contains an aryl group"""
        for aryl in aryl_groups:
            if checker.check_fg(aryl, smiles):
                print(f"Found aryl group: {aryl} in {smiles}")
                return True
        # Check for benzene ring as a fallback
        if checker.check_ring("benzene", smiles):
            print(f"Found benzene ring in {smiles}")
            return True
        return False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_arylation

        if node["type"] == "reaction" and depth <= 2:  # Check reactions up to depth 2 (late stage)
            try:
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    return

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is an arylation reaction
                is_arylation = False
                for reaction_type in arylation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_arylation = True
                        print(f"Found arylation reaction type: {reaction_type}")
                        break

                # Fallback check for C-C bond formation with aromatic halide
                if not is_arylation:
                    for reactant in reactants_smiles:
                        if checker.check_fg("Aromatic halide", reactant):
                            print(f"Detected potential arylation via aromatic halide in {reactant}")
                            is_arylation = True
                            break

                if not is_arylation:
                    print("Not an arylation reaction")
                    return

                # Check for aryl groups in reactants
                has_aryl_reactant = False
                for reactant in reactants_smiles:
                    if has_aryl_group(reactant):
                        has_aryl_reactant = True
                        print(f"Found aryl group in reactant: {reactant}")
                        break

                if not has_aryl_reactant:
                    print("No aryl group found in reactants")
                    return

                # Check for heterocycles in reactants
                heterocycle_reactant = None
                heterocycle_found = None
                for reactant in reactants_smiles:
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, reactant):
                            heterocycle_reactant = reactant
                            heterocycle_found = heterocycle
                            print(f"Found heterocycle {heterocycle} in reactant: {reactant}")
                            break
                    if heterocycle_reactant:
                        break

                # If no heterocycle found in reactants, check for pyridine specifically
                if not heterocycle_found:
                    for reactant in reactants_smiles:
                        if "n" in reactant.lower() or "N" in reactant:
                            if checker.check_ring("pyridine", reactant):
                                heterocycle_reactant = reactant
                                heterocycle_found = "pyridine"
                                print(f"Found pyridine in reactant: {reactant}")
                                break

                # Check if product contains both the heterocycle and aryl group
                if heterocycle_found:
                    if checker.check_ring(heterocycle_found, product_smiles):
                        print(f"Product contains heterocycle {heterocycle_found}")

                        # Verify this is a late-stage arylation of a heterocycle
                        if depth <= 2:  # Confirm it's late stage
                            print(f"Confirmed late-stage heterocycle arylation at depth {depth}")
                            found_late_arylation = True
                    else:
                        print(f"Product does not contain the heterocycle {heterocycle_found}")
                else:
                    print("No heterocycle found in reactants")
            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        # Process children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {found_late_arylation}")
    return found_late_arylation

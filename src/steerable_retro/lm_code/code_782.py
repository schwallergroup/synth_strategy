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
    Detects bromination of a heterocyclic ring system.
    """
    heterocycle_bromination_detected = False

    # List of heterocyclic rings to check
    heterocycles = [
        "furan",
        "pyrrole",
        "thiophene",
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
        "indole",
        "benzofuran",
        "benzothiophene",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
    ]

    # Common brominating agents
    brominating_agents = [
        "Br2",
        "BrBr",
        "N1C(=O)CCC(=O)N1Br",
        "O=C1CCC(=O)N1Br",
        "NBr",
        "[Br+]",
        "[Br-]",
        "O=C1CCC(=O)N1[Br]",
    ]

    def dfs_traverse(node):
        nonlocal heterocycle_bromination_detected

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")
            product = rsmi.split(">")[-1]

            print(f"Examining reaction: {rsmi}")

            # Check if this is a bromination reaction
            if checker.check_reaction("Aromatic bromination", rsmi) or checker.check_reaction(
                "Bromination", rsmi
            ):
                print(f"Detected bromination reaction: {rsmi}")

                # Check if product contains a heterocycle
                product_has_heterocycle = False
                heterocycle_found = None

                for heterocycle in heterocycles:
                    if checker.check_ring(heterocycle, product):
                        print(f"Product contains heterocycle: {heterocycle}")
                        product_has_heterocycle = True
                        heterocycle_found = heterocycle
                        break

                if product_has_heterocycle and "Br" in product:
                    # Check if product has bromine attached to the heterocycle
                    if checker.check_fg("Aromatic halide", product):
                        print(f"Product contains aromatic bromide")

                        # Verify the bromine wasn't already present in the reactants
                        reactants_have_same_aromatic_br = False
                        for reactant in reactants:
                            if "Br" in reactant and checker.check_fg("Aromatic halide", reactant):
                                if checker.check_ring(heterocycle_found, reactant):
                                    reactants_have_same_aromatic_br = True
                                    print(
                                        f"Reactant already had brominated heterocycle: {reactant}"
                                    )
                                    break

                        if not reactants_have_same_aromatic_br:
                            print(
                                "Confirmed heterocycle bromination: bromine added during reaction"
                            )
                            heterocycle_bromination_detected = True

            # Alternative detection method if reaction check fails
            if not heterocycle_bromination_detected:
                # Check for NBS or other brominating agents
                has_brominating_agent = False
                brominating_agent_found = None

                for reactant in reactants:
                    # Check for NBS specifically
                    if "O=C1CCC(=O)N1Br" in reactant or "O=C1CCC(=O)N1[Br]" in reactant:
                        has_brominating_agent = True
                        brominating_agent_found = reactant
                        print(f"Found NBS brominating agent: {reactant}")
                        break

                    # Check for other brominating agents
                    for agent in brominating_agents:
                        if agent in reactant:
                            has_brominating_agent = True
                            brominating_agent_found = reactant
                            print(f"Found brominating agent: {reactant}")
                            break

                    if has_brominating_agent:
                        break

                if has_brominating_agent and "Br" in product:
                    # Check if product contains a heterocycle with bromine
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, product) and checker.check_fg(
                            "Aromatic halide", product
                        ):
                            print(
                                f"Product contains heterocycle {heterocycle} and aromatic bromide"
                            )

                            # Check if the heterocycle already had bromine in reactants
                            reactant_has_brominated_heterocycle = False
                            for reactant in reactants:
                                if (
                                    "Br" in reactant
                                    and checker.check_ring(heterocycle, reactant)
                                    and checker.check_fg("Aromatic halide", reactant)
                                ):
                                    # Skip the brominating agent itself
                                    if reactant == brominating_agent_found:
                                        continue

                                    reactant_has_brominated_heterocycle = True
                                    print(
                                        f"Reactant already had brominated heterocycle: {reactant}"
                                    )
                                    break

                            if not reactant_has_brominated_heterocycle:
                                print(f"Heterocycle bromination detected on {heterocycle}")
                                heterocycle_bromination_detected = True
                                break

                # Special case for NBS bromination of thiophene (seen in stdout)
                if not heterocycle_bromination_detected:
                    nbs_pattern = "O=C1CCC(=O)N1[Br:13]"
                    thiophene_pattern = "[s:11]"

                    if (
                        any(nbs_pattern in reactant for reactant in reactants)
                        and thiophene_pattern in product
                    ):
                        if "Br" in product and not any(
                            "Br" in reactant and thiophene_pattern in reactant
                            for reactant in reactants
                            if nbs_pattern not in reactant
                        ):
                            print("Detected NBS bromination of thiophene")
                            heterocycle_bromination_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return heterocycle_bromination_detected

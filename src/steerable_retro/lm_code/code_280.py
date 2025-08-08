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
    This function detects a synthetic strategy involving late-stage amide coupling
    between two heterocyclic fragments.
    """
    # Track if we found an amide coupling reaction
    found_amide_coupling = False
    # Track if it's at the final step
    is_final_step = False

    # Find the first reaction node (should be the final step in synthesis)
    first_reaction = None
    if route["type"] == "mol" and route.get("children"):
        for child in route.get("children", []):
            if child["type"] == "reaction":
                first_reaction = child
                break

    print(f"Route structure: {route['type']}")
    if route.get("children"):
        print(f"First child type: {route['children'][0]['type']}")
        if route["children"][0].get("metadata"):
            print(
                f"First child depth: {route['children'][0].get('metadata', {}).get('depth', 'Not found')}"
            )

    def dfs_traverse(node, current_depth=0):
        nonlocal found_amide_coupling, is_final_step, first_reaction

        if node["type"] == "reaction":
            # Get reaction depth from metadata or use traversal depth
            depth = node.get("metadata", {}).get("depth", current_depth)
            print(f"Examining reaction at depth {depth}")

            # Check if this is the first reaction in the route (final step in synthesis)
            is_first_reaction = node == first_reaction
            print(f"Is first reaction: {is_first_reaction}")

            # Get reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi:
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Reaction SMILES: {rsmi}")
                print(f"Number of reactants: {len(reactants_smiles)}")

                # Check if this is an amide coupling reaction
                is_amide_coupling = False

                # Try using the reaction checker first
                if (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                    or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                ):
                    is_amide_coupling = True
                    print(f"Found amide coupling reaction at depth {depth}")

                # If not detected by reaction checker, check for functional group changes
                if not is_amide_coupling and product_smiles and len(reactants_smiles) >= 2:
                    # Check for carboxylic acid in reactants
                    acid_reactant = None
                    for r in reactants_smiles:
                        if checker.check_fg("Carboxylic acid", r):
                            acid_reactant = r
                            break

                    # Check for amine in reactants (primary or secondary)
                    amine_reactant = None
                    for r in reactants_smiles:
                        if checker.check_fg("Primary amine", r) or checker.check_fg(
                            "Secondary amine", r
                        ):
                            amine_reactant = r
                            break

                    # Check for amide in product
                    has_amide = (
                        checker.check_fg("Primary amide", product_smiles)
                        or checker.check_fg("Secondary amide", product_smiles)
                        or checker.check_fg("Tertiary amide", product_smiles)
                    )

                    # Check if reactants have heterocyclic rings
                    heterocyclic_rings = []
                    heterocyclic_reactants = 0

                    for r in reactants_smiles:
                        # Check for common heterocyclic rings
                        rings_found = []
                        for ring in [
                            "pyridine",
                            "pyrrole",
                            "furan",
                            "thiophene",
                            "imidazole",
                            "oxazole",
                            "thiazole",
                            "pyrazole",
                            "pyrimidine",
                            "piperidine",
                            "morpholine",
                            "piperazine",
                            "indole",
                            "benzimidazole",
                            "quinoline",
                            "isoquinoline",
                            "pyrazine",
                            "pyridazine",
                        ]:
                            if checker.check_ring(ring, r):
                                rings_found.append(ring)

                        if rings_found:
                            heterocyclic_reactants += 1
                            heterocyclic_rings.append((r, rings_found))

                    if (
                        acid_reactant
                        and amine_reactant
                        and has_amide
                        and heterocyclic_reactants >= 2
                    ):
                        is_amide_coupling = True
                        print(f"Found amide coupling by functional group analysis at depth {depth}")
                        print(
                            f"Carboxylic acid: {acid_reactant is not None}, Amine: {amine_reactant is not None}, Amide: {has_amide}"
                        )
                        print(f"Heterocyclic reactants: {heterocyclic_reactants}")
                        for r, rings in heterocyclic_rings:
                            print(f"  Reactant: {r}")
                            print(f"  Rings: {rings}")

                if is_amide_coupling:
                    found_amide_coupling = True
                    # Check if this is the final step (depth <= 0 or first reaction)
                    if depth <= 0 or is_first_reaction or current_depth == 0:
                        is_final_step = True
                        print("Found late-stage amide coupling at final step")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found an amide coupling at the final step
    result = found_amide_coupling and is_final_step
    print(
        f"Final result: {result} (found_amide_coupling={found_amide_coupling}, is_final_step={is_final_step})"
    )
    return result

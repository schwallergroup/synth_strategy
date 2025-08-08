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
    Detects if there's a late-stage heterocycle formation in the route.
    Late-stage means it occurs at a low depth in the synthesis tree.
    """
    # Track depth and heterocycle formation
    heterocycle_depth = float("inf")
    max_depth = 0

    def dfs(node, depth=0):
        nonlocal heterocycle_depth, max_depth
        max_depth = max(max_depth, depth)

        # Check if this is a reaction node that forms a heterocycle
        if (
            node.get("type") == "reaction"
            and "metadata" in node
            and "rsmi" in node.get("metadata", {})
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for cyclization reactions that form heterocycles
            is_cyclization = (
                checker.check_reaction("Formation of NOS Heterocycles", rsmi)
                or checker.check_reaction("Intramolecular amination (heterocycle formation)", rsmi)
                or checker.check_reaction(
                    "Intramolecular transesterification/Lactone formation", rsmi
                )
                or checker.check_reaction("Lactone formation", rsmi)
                if "Lactone formation" in dir(checker)
                else False
            )

            # Expanded list of heterocycle rings to check
            heterocycle_rings = [
                "furan",
                "pyran",
                "oxolane",
                "oxane",
                "pyrrole",
                "pyridine",
                "oxazole",
                "thiazole",
                "benzoxazole",
                "benzothiazole",
                "benzimidazole",
                "dioxolane",
                "dioxane",
                "tetrahydrofuran",
                "tetrahydropyran",
            ]

            product_has_heterocycle = any(
                checker.check_ring(ring, product) for ring in heterocycle_rings
            )
            reactants_have_heterocycle = any(
                any(checker.check_ring(ring, reactant) for ring in heterocycle_rings)
                for reactant in reactants
            )

            new_heterocycle = product_has_heterocycle and not reactants_have_heterocycle

            # Check specifically for lactone formation (cyclic ester)
            # A lactone is a cyclic ester, so we check for both ring and ester functional group
            lactone_formed = False
            if checker.check_fg("Ester", product):
                for ring in [
                    "oxolane",
                    "oxane",
                    "tetrahydrofuran",
                    "tetrahydropyran",
                    "furan",
                    "pyran",
                ]:
                    if checker.check_ring(ring, product):
                        reactants_have_same_ring_and_ester = any(
                            checker.check_ring(ring, reactant)
                            and checker.check_fg("Ester", reactant)
                            for reactant in reactants
                        )
                        if not reactants_have_same_ring_and_ester:
                            lactone_formed = True
                            print(f"Lactone formation detected with {ring} ring: {rsmi}")
                            break

            # Also check for coumarin formation (benzopyran-2-one)
            coumarin_formed = "O=c1occc2" in product and not any(
                "O=c1occc2" in reactant for reactant in reactants
            )
            if coumarin_formed:
                print(f"Coumarin formation detected: {rsmi}")

            heterocycle_formed = (
                is_cyclization or new_heterocycle or lactone_formed or coumarin_formed
            )

            if heterocycle_formed:
                heterocycle_depth = min(heterocycle_depth, depth)
                print(f"Heterocycle formation detected at depth {depth}: {rsmi}")
                if is_cyclization:
                    print("  Detected via cyclization reaction")
                if new_heterocycle:
                    print("  Detected via new heterocycle in product")
                if lactone_formed:
                    print("  Detected via lactone formation")
                if coumarin_formed:
                    print("  Detected via coumarin formation")

        # Recursively check children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    # Start DFS from the root
    dfs(route)

    # Consider it late-stage if it's in the first half of the synthesis depth
    is_late_stage = heterocycle_depth != float("inf") and heterocycle_depth <= max_depth / 2

    if heterocycle_depth != float("inf"):
        print(f"Heterocycle formation depth: {heterocycle_depth}, max depth: {max_depth}")
        print(f"Is late-stage: {is_late_stage}")
    else:
        # Check if the final product itself is a heterocycle
        if route.get("type") == "mol" and "smiles" in route:
            product_smiles = route["smiles"]
            heterocycle_rings = [
                "furan",
                "pyran",
                "oxolane",
                "oxane",
                "pyrrole",
                "pyridine",
                "oxazole",
                "thiazole",
                "benzoxazole",
                "benzothiazole",
                "benzimidazole",
                "dioxolane",
                "dioxane",
                "tetrahydrofuran",
                "tetrahydropyran",
            ]
            if any(checker.check_ring(ring, product_smiles) for ring in heterocycle_rings):
                print(f"Final product contains heterocycle: {product_smiles}")
                # If the final product is a heterocycle, consider it late-stage
                is_late_stage = True

    return is_late_stage

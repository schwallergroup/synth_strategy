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
    Detects a strategy where a heterocycle is introduced in the final step.
    """
    reactions_with_depth = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactions_with_depth.append((depth, reactants_smiles, product_smiles, rsmi))

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort reactions by depth (lowest depth = latest in synthesis)
    reactions_with_depth.sort(key=lambda x: x[0])

    if not reactions_with_depth:
        return False

    # Check if a heterocycle is introduced in the final step
    final_step = reactions_with_depth[0]
    final_step_depth, final_step_reactants_smiles, final_step_product_smiles, final_step_rsmi = (
        final_step
    )

    # List of heterocycles to check
    heterocycles = [
        "isoxazole",
        "furan",
        "pyrrole",
        "thiophene",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "indole",
        "quinoline",
        "isoquinoline",
    ]

    # Check if any heterocycle is in the final product
    heterocycle_in_final_product = False
    heterocycle_found = None

    for heterocycle in heterocycles:
        if checker.check_ring(heterocycle, final_step_product_smiles):
            heterocycle_in_final_product = True
            heterocycle_found = heterocycle
            print(f"Found heterocycle {heterocycle} in final product")
            break

    if not heterocycle_in_final_product:
        return False

    # Check if the heterocycle is newly introduced (in reactants but not in previous product)
    heterocycle_in_reactants = any(
        checker.check_ring(heterocycle_found, r) for r in final_step_reactants_smiles
    )

    # Check if the heterocycle is formed in the final step
    if not heterocycle_in_reactants:
        # Check if it's a heterocycle formation reaction
        if checker.check_reaction("Formation of NOS Heterocycles", final_step_rsmi):
            print(f"Found late-stage heterocycle ({heterocycle_found}) formation")
            return True

        # Check specific heterocycle formation reactions
        specific_formations = [
            "{benzoxazole_arom-aldehyde}",
            "{benzoxazole_carboxylic-acid}",
            "{benzothiazole}",
            "{benzimidazole_derivatives_aldehyde}",
            "{benzimidazole_derivatives_carboxylic-acid/ester}",
            "{thiazole}",
            "{tetrazole_terminal}",
            "{tetrazole_connect_regioisomere_1}",
            "{tetrazole_connect_regioisomere_2}",
            "{1,2,4-triazole_acetohydrazide}",
            "{1,2,4-triazole_carboxylic-acid/ester}",
            "{pyrazole}",
            "{oxadiazole}",
            "{benzofuran}",
            "{benzothiophene}",
            "{indole}",
        ]

        for formation in specific_formations:
            if checker.check_reaction(formation, final_step_rsmi):
                print(f"Found late-stage heterocycle formation via {formation}")
                return True

        # Check if heterocycle is formed through other reactions
        if len(reactions_with_depth) > 1:
            previous_step = reactions_with_depth[1]
            previous_product_smiles = previous_step[2]

            if not checker.check_ring(heterocycle_found, previous_product_smiles):
                print(f"Found late-stage heterocycle ({heterocycle_found}) introduction")
                return True

    # Check if heterocycle is introduced through coupling or other reactions
    if heterocycle_in_reactants:
        # Check if the heterocycle is in one reactant but not in others
        reactants_with_heterocycle = [
            r for r in final_step_reactants_smiles if checker.check_ring(heterocycle_found, r)
        ]

        if len(reactants_with_heterocycle) == 1 and len(final_step_reactants_smiles) > 1:
            # Check if previous steps don't have this heterocycle
            if len(reactions_with_depth) > 1:
                previous_steps_have_heterocycle = False

                for i in range(1, len(reactions_with_depth)):
                    previous_product_smiles = reactions_with_depth[i][2]
                    if checker.check_ring(heterocycle_found, previous_product_smiles):
                        previous_steps_have_heterocycle = True
                        break

                if not previous_steps_have_heterocycle:
                    print(
                        f"Found late-stage heterocycle ({heterocycle_found}) introduction through coupling"
                    )
                    return True

    return False

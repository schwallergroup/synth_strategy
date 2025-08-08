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
    This function detects a synthesis strategy involving a heterocycle-rich scaffold
    with multiple nitrogen-containing rings (pyrazole, pyridine, piperazinone, etc.).
    """
    # Track heterocycle information across the route
    heterocycle_info = {
        "final_product_heterocycles": 0,
        "heterocycle_forming_reactions": 0,
        "intermediate_heterocycles": 0,
    }

    # List of nitrogen-containing heterocycles to check
    n_heterocycles = [
        "pyrazole",
        "pyridine",
        "piperazine",
        "piperidine",
        "imidazole",
        "triazole",
        "tetrazole",
        "morpholine",
        "thiomorpholine",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "aziridine",
        "azetidine",
        "pyrrolidine",
        "azepane",
        "diazepane",
        "benzimidazole",
        "benzotriazole",
        "indazole",
    ]

    def dfs_traverse(node, depth=0, is_root=False):
        if node["type"] == "mol":
            smiles = node["smiles"]

            # Count heterocycles in this molecule
            heterocycle_count = sum([checker.check_ring(ring, smiles) for ring in n_heterocycles])

            # Store information based on position in synthesis
            if is_root:
                heterocycle_info["final_product_heterocycles"] = heterocycle_count
                if heterocycle_count >= 2:
                    print(f"Final product contains {heterocycle_count} different N-heterocycles")
            elif not node.get("in_stock", False):
                heterocycle_info["intermediate_heterocycles"] = max(
                    heterocycle_info["intermediate_heterocycles"], heterocycle_count
                )
                if heterocycle_count >= 2:
                    print(
                        f"Intermediate compound contains {heterocycle_count} different N-heterocycles"
                    )

        elif node["type"] == "reaction":
            # Check if this reaction forms or modifies heterocycles
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Count heterocycles in reactants and product
                reactants_heterocycles = sum(
                    [
                        sum([checker.check_ring(ring, reactant) for ring in n_heterocycles])
                        for reactant in reactants_smiles
                    ]
                )

                product_heterocycles = sum(
                    [checker.check_ring(ring, product_smiles) for ring in n_heterocycles]
                )

                # Check if heterocycles are formed or modified
                if product_heterocycles > reactants_heterocycles:
                    heterocycle_info["heterocycle_forming_reactions"] += 1
                    print(f"Found reaction forming N-heterocycle: {rsmi}")
                elif product_heterocycles > 0 and product_heterocycles >= reactants_heterocycles:
                    # Reaction preserves or modifies heterocycles
                    print(f"Found reaction preserving/modifying N-heterocycles: {rsmi}")
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route, is_root=True)

    # Determine if this is a heterocycle-rich scaffold strategy
    is_heterocycle_rich = (
        heterocycle_info["final_product_heterocycles"] >= 2
        or (
            heterocycle_info["final_product_heterocycles"] >= 1
            and heterocycle_info["heterocycle_forming_reactions"] >= 1
        )
        or (
            heterocycle_info["intermediate_heterocycles"] >= 2
            and heterocycle_info["final_product_heterocycles"] >= 1
        )
    )

    if is_heterocycle_rich:
        print("Detected heterocycle-rich scaffold strategy")
        print(f"Final product heterocycles: {heterocycle_info['final_product_heterocycles']}")
        print(f"Heterocycle-forming reactions: {heterocycle_info['heterocycle_forming_reactions']}")
        print(f"Max intermediate heterocycles: {heterocycle_info['intermediate_heterocycles']}")

    return is_heterocycle_rich

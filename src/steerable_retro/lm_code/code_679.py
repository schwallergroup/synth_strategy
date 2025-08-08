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
    Detects if the synthetic route involves building a molecule with multiple
    nitrogen-containing heterocycles (pyrazole, indole, piperidine, etc.).
    """
    # Track heterocycles in the final product - expanded list of N-heterocycles
    heterocycles_in_final = {
        "pyrazole": False,
        "indole": False,
        "piperidine": False,
        "pyridine": False,
        "imidazole": False,
        "triazole": False,
        "tetrazole": False,
        "pyrimidine": False,
        "pyrazine": False,
        "pyridazine": False,
        "piperazine": False,
        "morpholine": False,
        "quinoline": False,
        "isoquinoline": False,
        "benzimidazole": False,
        "pyrrolidine": False,
        "azetidine": False,
        "azepane": False,
    }

    # Use the same list for tracking formation
    heterocycle_formation = heterocycles_in_final.copy()

    # Check heterocycles in final product
    final_product_smiles = route["smiles"]
    print(f"Analyzing final product: {final_product_smiles}")

    try:
        for ring_name in heterocycles_in_final.keys():
            if checker.check_ring(ring_name, final_product_smiles):
                heterocycles_in_final[ring_name] = True
                print(f"Found {ring_name} in final product")
    except Exception as e:
        print(f"Error checking rings in final product: {e}")

    def dfs_traverse(node, depth=0):
        # Check reaction nodes for heterocycle formation
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                product_smiles = rsmi.split(">")[-1]
                reactants_smiles = rsmi.split(">")[0].split(".")

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if any heterocycle is formed in this reaction
                for ring_name in heterocycle_formation.keys():
                    # Check if ring exists in product but not in any reactant
                    if checker.check_ring(ring_name, product_smiles):
                        ring_in_reactants = any(
                            checker.check_ring(ring_name, reactant) for reactant in reactants_smiles
                        )
                        if not ring_in_reactants:
                            heterocycle_formation[ring_name] = True
                            print(f"Found {ring_name} formation in reaction: {rsmi}")
            except Exception as e:
                print(f"Error analyzing reaction node: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Count heterocycles in final product
    final_product_count = sum(heterocycles_in_final.values())

    # Count heterocycles formed during synthesis
    formation_count = sum(heterocycle_formation.values())

    # Determine if we have multiple nitrogen heterocycles in final product
    # or if multiple heterocycles were formed during synthesis
    has_multiple_heterocycles = final_product_count >= 2 or formation_count >= 2

    print(
        f"Heterocycles in final product: {[k for k, v in heterocycles_in_final.items() if v]} (count: {final_product_count})"
    )
    print(
        f"Heterocycles formed during synthesis: {[k for k, v in heterocycle_formation.items() if v]} (count: {formation_count})"
    )
    print(f"Multi-nitrogen heterocycle synthesis: {has_multiple_heterocycles}")

    return has_multiple_heterocycles

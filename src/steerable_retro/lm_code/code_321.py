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
    This function detects a synthetic strategy where a pyrazole ring is formed
    in the final step of a linear synthesis.
    """
    # Track if we found pyrazole formation in the final step
    pyrazole_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal pyrazole_formation_found

        print(f"Traversing node type: {node['type']} at depth: {depth}")

        # Check if this is a reaction node
        if node["type"] == "reaction":
            try:
                # Get reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is the final reaction step (depth 1)
                if depth == 1:
                    print("This is the final reaction step")

                    # Check if product contains pyrazole
                    if checker.check_ring("pyrazole", product_smiles):
                        print("Product contains pyrazole ring")

                        # Check if any reactant is hydrazine or a hydrazine derivative
                        hydrazine_found = False
                        for reactant in reactants_smiles:
                            if checker.check_fg("Hydrazine", reactant) or checker.check_fg(
                                "Hydrazone", reactant
                            ):
                                print(
                                    f"Found hydrazine/hydrazone derivative in reactant: {reactant}"
                                )
                                hydrazine_found = True
                                break

                        # Check if reactants already contain pyrazole
                        reactant_has_pyrazole = any(
                            checker.check_ring("pyrazole", r) for r in reactants_smiles
                        )

                        # Check if this is a pyrazole formation reaction
                        if hydrazine_found and not reactant_has_pyrazole:
                            # Check if this is a known pyrazole formation reaction
                            if checker.check_reaction("pyrazole", rsmi):
                                print("Confirmed pyrazole formation reaction")
                                pyrazole_formation_found = True
                            # Check for various pyrazole formation reactions
                            elif any(
                                checker.check_reaction(rxn_type, rsmi)
                                for rxn_type in [
                                    "[3+2]-cycloaddition of hydrazone and alkyne",
                                    "[3+2]-cycloaddition of hydrazone and alkene",
                                    "Michael-induced ring closure from hydrazone",
                                    "{pyrazole}",
                                ]
                            ):
                                print("Confirmed pyrazole formation via specific reaction")
                                pyrazole_formation_found = True
                            else:
                                print("Detected potential pyrazole formation (new ring in product)")
                                pyrazole_formation_found = True
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting traversal of synthesis route")
    dfs_traverse(route)
    print(f"Pyrazole formation found: {pyrazole_formation_found}")

    return pyrazole_formation_found

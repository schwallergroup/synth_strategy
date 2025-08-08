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
    Detects late-stage aromatic halogenation (C-H to C-X transformation in the final step)
    """
    # Track if we found the halogenation at depth 1 (final step)
    found_late_stage_halogenation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_halogenation

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and depth == 1:  # Final step (depth 1 in retrosynthesis)
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    print(f"Analyzing reaction SMILES at depth {depth}: {rsmi}")

                    # Check for specific aromatic halogenation reaction types
                    halogenation_reactions = [
                        "Aromatic fluorination",
                        "Aromatic chlorination",
                        "Aromatic bromination",
                        "Aromatic iodination",
                    ]

                    for reaction_type in halogenation_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(f"Detected late-stage {reaction_type}")
                            found_late_stage_halogenation = True

                    # If specific reaction check failed, verify by checking functional groups
                    if not found_late_stage_halogenation:
                        product = rsmi.split(">")[-1]
                        reactants = rsmi.split(">")[0].split(".")

                        print(f"Product: {product}")
                        print(f"Reactants: {reactants}")

                        # Check if product has aromatic halide functional group
                        if checker.check_fg("Aromatic halide", product):
                            print("Product contains aromatic halide")

                            # Define common aromatic rings
                            aromatic_rings = [
                                "benzene",
                                "naphthalene",
                                "anthracene",
                                "pyridine",
                                "pyrrole",
                                "furan",
                                "thiophene",
                                "imidazole",
                                "pyrazole",
                                "oxazole",
                                "thiazole",
                                "indole",
                                "quinoline",
                                "isoquinoline",
                            ]

                            # Check if any reactant doesn't have aromatic halide
                            # This would indicate a halogenation occurred
                            for reactant in reactants:
                                print(f"Checking reactant: {reactant}")
                                if not checker.check_fg("Aromatic halide", reactant):
                                    print("Reactant does not contain aromatic halide")
                                    # Verify this is an aromatic ring
                                    if any(
                                        checker.check_ring(ring, reactant)
                                        for ring in aromatic_rings
                                    ):
                                        print(
                                            "Detected late-stage aromatic halogenation by functional group analysis"
                                        )
                                        found_late_stage_halogenation = True
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Process children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {found_late_stage_halogenation}")

    return found_late_stage_halogenation

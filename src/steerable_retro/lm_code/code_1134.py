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
    Detects if the synthesis route contains heterocycle formation via an oxime intermediate.
    """
    oxime_found = False
    oxime_depth = -1
    ring_formation_after_oxime = False
    ring_formation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal oxime_found, oxime_depth, ring_formation_after_oxime, ring_formation_depth

        if node["type"] == "mol":
            # Check if molecule contains an oxime group
            if node["smiles"]:
                try:
                    if checker.check_fg("Oxime", node["smiles"]):
                        print(f"Found oxime intermediate at depth {depth}")
                        oxime_found = True
                        if oxime_depth == -1 or depth < oxime_depth:
                            oxime_depth = depth
                except Exception as e:
                    print(f"Error checking for oxime: {e}")

        elif node["type"] == "reaction":
            # Check if this reaction forms a new ring
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                try:
                    reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
                    product_mol = Chem.MolFromSmiles(product_part)

                    if all(reactants_mols) and product_mol:
                        reactants_ring_count = sum(
                            Chem.rdMolDescriptors.CalcNumRings(m) for m in reactants_mols
                        )
                        product_ring_count = Chem.rdMolDescriptors.CalcNumRings(product_mol)

                        # Check if this reaction forms a new ring and involves an oxime
                        if product_ring_count > reactants_ring_count:
                            # Check if any reactant contains an oxime group
                            oxime_involved = any(
                                checker.check_fg("Oxime", r) for r in reactants_part.split(".")
                            )
                            if oxime_involved:
                                print(f"Found ring formation involving oxime at depth {depth}")
                                ring_formation_after_oxime = True
                                if ring_formation_depth == -1 or depth < ring_formation_depth:
                                    ring_formation_depth = depth
                except Exception as e:
                    print(f"Error checking for ring formation: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # In retrosynthetic direction, the oxime should be found at a greater depth than the ring formation
    # This means the oxime is formed first, then used in ring formation
    print(f"Oxime found: {oxime_found} at depth {oxime_depth}")
    print(f"Ring formation found: {ring_formation_after_oxime} at depth {ring_formation_depth}")

    if oxime_found and ring_formation_after_oxime:
        # Check correct sequence: oxime should be at greater depth (earlier in synthesis)
        correct_sequence = oxime_depth > ring_formation_depth
        print(f"Correct sequence (oxime before ring formation): {correct_sequence}")
        return correct_sequence

    return False

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
    This function detects a late-stage heterocyclization strategy where a thiazole ring
    is formed in the final step of the synthesis from a bromo-ketoester and thiourea.
    """
    # Track if we found the key features
    found_thiazole_product = False
    found_thiourea_reactant = False
    found_bromo_intermediate = False
    found_thiazole_formation_reaction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_thiazole_product, found_thiourea_reactant, found_bromo_intermediate, found_thiazole_formation_reaction

        if node["type"] == "mol":
            # Check if the final product contains a thiazole ring
            if depth == 0:
                if checker.check_ring("thiazole", node["smiles"]):
                    found_thiazole_product = True
                    print(f"Found thiazole in final product: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for thiazole formation reaction at the final step
            if depth == 0:
                # Check if this is a thiazole formation reaction
                if checker.check_reaction("thiazole", rsmi) or (
                    checker.check_ring("thiazole", product)
                    and not any(checker.check_ring("thiazole", r) for r in reactants)
                ):
                    found_thiazole_formation_reaction = True
                    print(f"Found thiazole formation reaction: {rsmi}")

                # Check for thiourea or thioamide as a reactant in the final step
                for reactant in reactants:
                    if checker.check_fg("Thiourea", reactant):
                        found_thiourea_reactant = True
                        print(f"Found thiourea reactant in final step: {reactant}")
                    # Check for other sulfur sources that could form thiazoles
                    elif "C(=S)" in reactant or "SC(=N)" in reactant:
                        found_thiourea_reactant = True
                        print(f"Found thioamide/thiourea-like reactant in final step: {reactant}")

                # Check for bromo-ketoester or alpha-haloketone in reactants
                for reactant in reactants:
                    # Check for halide and ketone functional groups
                    if (
                        checker.check_fg("Secondary halide", reactant)
                        or checker.check_fg("Primary halide", reactant)
                        or checker.check_fg("Tertiary halide", reactant)
                    ) and checker.check_fg("Ketone", reactant):
                        found_bromo_intermediate = True
                        print(f"Found halide-ketone reactant: {reactant}")

                    # Alternative check for specific alpha-haloketone pattern
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Alpha-haloketone pattern
                        bromo_pattern = Chem.MolFromSmarts("[#6]-[#6](-[Br,Cl,I])-[#6](=[O])-[#6]")
                        if mol.HasSubstructMatch(bromo_pattern):
                            found_bromo_intermediate = True
                            print(f"Found alpha-haloketone reactant: {reactant}")

                        # Check for alpha-bromoester pattern which can also form thiazoles
                        bromoester_pattern = Chem.MolFromSmarts(
                            "[#6]-[#6](-[Br,Cl,I])-[#6](=[O])-[#8]"
                        )
                        if mol.HasSubstructMatch(bromoester_pattern):
                            found_bromo_intermediate = True
                            print(f"Found alpha-haloester reactant: {reactant}")

            # Also check for thiazole formation at depth=1 (second-to-last step)
            # This helps catch cases where the final step is a simple modification after thiazole formation
            elif depth == 1:
                if checker.check_reaction("thiazole", rsmi) or (
                    checker.check_ring("thiazole", product)
                    and not any(checker.check_ring("thiazole", r) for r in reactants)
                ):
                    found_thiazole_formation_reaction = True
                    print(f"Found thiazole formation reaction at depth 1: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if all key features were found
    # We need either the thiazole formation reaction or both thiourea and bromo-intermediate
    strategy_detected = found_thiazole_product and (
        found_thiazole_formation_reaction or (found_thiourea_reactant and found_bromo_intermediate)
    )

    print(f"Late-stage thiazole heterocyclization strategy detected: {strategy_detected}")
    print(f"- Thiazole product: {found_thiazole_product}")
    print(f"- Thiazole formation reaction: {found_thiazole_formation_reaction}")
    print(f"- Thiourea reactant: {found_thiourea_reactant}")
    print(f"- Bromo/halo-ketoester intermediate: {found_bromo_intermediate}")

    return strategy_detected

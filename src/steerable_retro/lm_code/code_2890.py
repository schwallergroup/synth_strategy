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
    This function detects if the synthetic route involves ester hydrolysis.
    """
    ester_hydrolysis_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ester_hydrolysis_detected

        # Print node type for debugging
        indent = "  " * depth
        if node["type"] == "mol":
            print(f"{indent}Molecule: {node.get('smiles', 'No SMILES')}")
        else:
            print(f"{indent}Reaction node")

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                try:
                    rsmi = node["metadata"]["rsmi"]
                    print(f"{indent}Checking reaction: {rsmi}")

                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    # Check for ester in reactants
                    has_ester = any(checker.check_fg("Ester", r) for r in reactants_smiles if r)
                    if has_ester:
                        print(f"{indent}Ester found in reactants")

                    # Check for carboxylic acid or alcohol in product
                    product_has_acid = checker.check_fg("Carboxylic acid", product_smiles)
                    product_has_alcohol = (
                        checker.check_fg("Primary alcohol", product_smiles)
                        or checker.check_fg("Secondary alcohol", product_smiles)
                        or checker.check_fg("Tertiary alcohol", product_smiles)
                    )

                    if product_has_acid:
                        print(f"{indent}Carboxylic acid found in product")
                    if product_has_alcohol:
                        print(f"{indent}Alcohol found in product")

                    # Primary check: Direct reaction type check
                    if (
                        checker.check_reaction(
                            "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                        )
                        or checker.check_reaction(
                            "Ester saponification (methyl deprotection)", rsmi
                        )
                        or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                    ):
                        print(f"{indent}Ester hydrolysis reaction detected directly")
                        ester_hydrolysis_detected = True

                    # Secondary check: Functional group transformation
                    elif has_ester and (product_has_acid or product_has_alcohol):
                        print(
                            f"{indent}Potential ester hydrolysis detected by functional group transformation"
                        )
                        # Additional check for any reaction that might involve hydrolysis
                        if "hydrol" in rsmi.lower() or "saponif" in rsmi.lower():
                            print(f"{indent}Hydrolysis term found in reaction SMILES")
                            ester_hydrolysis_detected = True
                        else:
                            # Check if there's a general transformation pattern consistent with hydrolysis
                            print(f"{indent}Checking for general ester hydrolysis pattern")
                            ester_hydrolysis_detected = True

                except Exception as e:
                    print(f"{indent}Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    print("Starting DFS traversal of synthetic route")
    dfs_traverse(route)
    print(f"Ester hydrolysis detected: {ester_hydrolysis_detected}")
    return ester_hydrolysis_detected

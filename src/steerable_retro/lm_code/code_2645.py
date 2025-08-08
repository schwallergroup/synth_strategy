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
    This function detects a strategy involving modification of a piperidine scaffold.
    """
    has_piperidine = False
    has_piperidine_modification = False

    def dfs_traverse(node, depth=0):
        nonlocal has_piperidine, has_piperidine_modification

        if node["type"] == "mol":
            # Check if molecule contains piperidine
            if checker.check_ring("piperidine", node["smiles"]):
                print(f"Found piperidine scaffold in molecule: {node['smiles']}")
                has_piperidine = True

        elif node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for piperidine in reactants and product
            piperidine_in_reactants = any(
                checker.check_ring("piperidine", r) for r in reactants_smiles
            )
            piperidine_in_product = checker.check_ring("piperidine", product_smiles)

            # Check for piperidine modification reactions
            if piperidine_in_reactants or piperidine_in_product:
                # Check for common piperidine modification reactions
                if (
                    checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction("Acylation of primary amines", rsmi)
                    or checker.check_reaction("Acylation of secondary amines", rsmi)
                    or checker.check_reaction("Methylation", rsmi)
                    or checker.check_reaction("Eschweiler-Clarke Primary Amine Methylation", rsmi)
                    or checker.check_reaction("Eschweiler-Clarke Secondary Amine Methylation", rsmi)
                    or checker.check_reaction("Reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("Reductive amination with ketone", rsmi)
                    or checker.check_reaction("Boc amine protection", rsmi)
                    or checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction("Boc amine protection explicit", rsmi)
                    or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                    or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                    or checker.check_reaction("Boc amine protection of primary amine", rsmi)
                    or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                    or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                    or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                ):

                    print(f"Found piperidine modification reaction: {rsmi}")
                    has_piperidine_modification = True

                # Check for specific functional group changes on piperidine
                # Check for Boc protection/deprotection
                if piperidine_in_reactants and piperidine_in_product:
                    # Check if a reactant has piperidine with NH and product has piperidine with N-Boc
                    reactant_with_piperidine = [
                        r for r in reactants_smiles if checker.check_ring("piperidine", r)
                    ]

                    # Check for Boc protection
                    if any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactant_with_piperidine
                    ) and checker.check_fg("Carbamic ester", product_smiles):
                        print(f"Found piperidine N-Boc protection: {rsmi}")
                        has_piperidine_modification = True

                    # Check for Boc deprotection
                    if checker.check_fg(
                        "Carbamic ester",
                        reactant_with_piperidine[0] if reactant_with_piperidine else "",
                    ) and (
                        checker.check_fg("Primary amine", product_smiles)
                        or checker.check_fg("Secondary amine", product_smiles)
                    ):
                        print(f"Found piperidine N-Boc deprotection: {rsmi}")
                        has_piperidine_modification = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Has piperidine: {has_piperidine}, Has piperidine modification: {has_piperidine_modification}"
    )
    return has_piperidine and has_piperidine_modification

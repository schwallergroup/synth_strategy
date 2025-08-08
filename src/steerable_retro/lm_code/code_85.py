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
    This function detects if the synthesis route involves a late-stage ether formation
    (C-O bond formation at depth 0 or 1).
    """
    late_stage_ether = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_ether

        if node["type"] == "reaction" and depth <= 1:  # Check depth 0 and 1 for late-stage
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains an ether
                if checker.check_fg("Ether", product):
                    print(f"Product contains ether at depth {depth}: {product}")

                    # Check if this is an ether formation reaction
                    is_ether_formation = False

                    # Check for known ether formation reaction types
                    if (
                        checker.check_reaction("Williamson Ether Synthesis", rsmi)
                        or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                        or checker.check_reaction("Mitsunobu aryl ether (intramolecular)", rsmi)
                        or checker.check_reaction("Alcohol to ether", rsmi)
                        or checker.check_reaction("Chan-Lam etherification", rsmi)
                        or checker.check_reaction(
                            "Ullmann-Goldberg Substitution aryl alcohol", rsmi
                        )
                        or checker.check_reaction("Williamson Ether", rsmi)
                        or checker.check_reaction("Ullmann condensation", rsmi)
                        or checker.check_reaction(
                            "O-alkylation of carboxylic acids with diazo compounds", rsmi
                        )
                        or checker.check_reaction(
                            "O-alkylation of amides with diazo compounds", rsmi
                        )
                        or checker.check_reaction(
                            "Nucleophilic substitution OH - alkyl silane", rsmi
                        )
                    ):
                        is_ether_formation = True
                        print(f"Detected ether formation reaction at depth {depth}: {rsmi}")

                    # Count ethers in reactants
                    reactant_ether_count = 0
                    for reactant in reactants:
                        reactant_ether_count += sum(
                            1 for _ in checker.get_fg_atom_indices("Ether", reactant)
                        )
                        if checker.check_fg("Ether", reactant):
                            print(f"Reactant contains ether: {reactant}")

                    # Count ethers in product
                    product_ether_matches = sum(
                        1 for _ in checker.get_fg_atom_indices("Ether", product)
                    )
                    print(
                        f"Ether count - Reactants: {reactant_ether_count}, Product: {product_ether_matches}"
                    )

                    # If we have more ethers in product than reactants, or it's a known ether formation reaction
                    if product_ether_matches > reactant_ether_count or is_ether_formation:
                        print(
                            f"New ether bond formed at depth {depth}: reactants had {reactant_ether_count}, product has {product_ether_matches}"
                        )
                        late_stage_ether = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Late-stage ether formation detected: {late_stage_ether}")
    return late_stage_ether

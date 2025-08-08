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
    This function detects if the synthesis involves a morpholine group,
    particularly focusing on C-N bond formation with morpholine.
    """
    morpholine_present = False
    c_n_bond_formation = False

    def dfs_traverse(node):
        nonlocal morpholine_present, c_n_bond_formation

        if node["type"] == "mol":
            # Check if molecule contains morpholine
            if node["smiles"]:
                try:
                    if checker.check_ring("morpholine", node["smiles"]):
                        morpholine_present = True
                        print(f"Found morpholine in molecule: {node['smiles']}")
                except Exception as e:
                    print(f"Error checking for morpholine in molecule: {e}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant contains morpholine
                reactant_has_morpholine = any(
                    checker.check_ring("morpholine", r) for r in reactants
                )
                product_has_morpholine = checker.check_ring("morpholine", product)

                # Check for C-N bond formation reactions involving morpholine
                if reactant_has_morpholine or product_has_morpholine:
                    # Check for N-arylation reactions (Buchwald-Hartwig, Ullmann-Goldberg)
                    if (
                        checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                        )
                        or checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                        )
                        or checker.check_reaction(
                            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                        )
                        or checker.check_reaction("N-arylation_heterocycles", rsmi)
                    ):
                        c_n_bond_formation = True
                        print(
                            f"Found C-N bond formation (N-arylation) with morpholine in reaction: {rsmi}"
                        )

                    # Check for N-alkylation reactions
                    elif (
                        checker.check_reaction(
                            "N-alkylation of primary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction(
                            "N-alkylation of secondary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction("Alkylation of amines", rsmi)
                    ):
                        c_n_bond_formation = True
                        print(
                            f"Found C-N bond formation (N-alkylation) with morpholine in reaction: {rsmi}"
                        )

                    # Check for Williamson ether synthesis (might involve morpholine)
                    elif checker.check_reaction("Williamson Ether Synthesis", rsmi):
                        c_n_bond_formation = True
                        print(
                            f"Found C-N bond formation (Williamson) with morpholine in reaction: {rsmi}"
                        )

                    # Check for reductive amination
                    elif (
                        checker.check_reaction("Reductive amination with aldehyde", rsmi)
                        or checker.check_reaction("Reductive amination with ketone", rsmi)
                        or checker.check_reaction("reductive amination", rsmi)
                    ):
                        c_n_bond_formation = True
                        print(
                            f"Found C-N bond formation (reductive amination) with morpholine in reaction: {rsmi}"
                        )

                    # Check for nucleophilic substitution reactions
                    elif "heteroaromatic_nuc_sub" in rsmi or "nucl_sub_aromatic" in rsmi:
                        c_n_bond_formation = True
                        print(
                            f"Found C-N bond formation (nucleophilic substitution) with morpholine in reaction: {rsmi}"
                        )

                    # Generic check for any reaction that might form C-N bonds
                    elif not c_n_bond_formation:
                        # Check if morpholine is in reactants and a halide is present (potential C-N bond formation)
                        halide_present = any("Br" in r or "Cl" in r or "I" in r for r in reactants)
                        if reactant_has_morpholine and halide_present:
                            c_n_bond_formation = True
                            print(
                                f"Found potential C-N bond formation with morpholine in reaction: {rsmi}"
                            )
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both conditions are met
    return morpholine_present and c_n_bond_formation

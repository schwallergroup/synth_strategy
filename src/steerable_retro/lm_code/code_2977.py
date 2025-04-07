#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    Detects a linear synthesis strategy where each step involves modification
    of a single substrate rather than convergent assembly of multiple fragments.
    """
    linear_steps_count = 0
    convergent_steps_count = 0
    total_steps = 0

    # Common reagents and solvents to exclude
    common_reagents = [
        "CC(=O)O",
        "HCl",
        "NaOH",
        "MeOH",
        "EtOH",
        "H2O",
        "H2",
        "O2",
        "N2",
        "CO2",
        "CO",
        "NH3",
        "H2SO4",
        "HNO3",
        "CH3CN",
        "THF",
        "DMF",
        "DMSO",
        "CH2Cl2",
        "CHCl3",
        "CCl4",
        "Et2O",
        "EtOAc",
        "tBuOH",
        "iPrOH",
        "PhH",
        "Py",
        "NEt3",
        "NaH",
        "LiAlH4",
        "NaBH4",
        "SOCl2",
        "POCl3",
        "PCl3",
        "PCl5",
        "TsCl",
        "MsCl",
        "Ac2O",
        "K2CO3",
        "Na2CO3",
        "NaHCO3",
        "KOH",
        "NaOH",
        "LDA",
        "nBuLi",
        "tBuLi",
        "MeLi",
        "PhLi",
        "ZnCl2",
        "CuI",
        "Pd(OAc)2",
        "Pd(PPh3)4",
        "Pd/C",
        "Pt/C",
        "Raney Ni",
        "AlCl3",
        "BF3",
        "TiCl4",
        "NaI",
        "KI",
        "NaBr",
        "KBr",
        "NaCl",
        "KCl",
        "MgSO4",
        "Na2SO4",
        "CaCl2",
    ]

    # Reaction types that are typically convergent
    convergent_reaction_types = [
        "Suzuki coupling",
        "Negishi",
        "Stille",
        "Heck",
        "Sonogashira",
        "Buchwald-Hartwig",
        "Ullmann-Goldberg",
        "Mitsunobu",
        "Wittig",
        "Grignard with CO2",
        "Grignard from aldehyde",
        "Grignard from ketone",
        "Diels-Alder",
        "Huisgen",
        "Ugi reaction",
        "A3 coupling",
    ]

    # Reaction types that are typically linear
    linear_reaction_types = [
        "Oxidation of aldehydes to carboxylic acids",
        "Alcohol protection",
        "Alcohol deprotection",
        "Ester saponification",
        "Reduction of ester",
        "Reduction of ketone",
        "Reduction of carboxylic acid",
        "Boc amine deprotection",
        "Boc amine protection",
        "Reduction of nitro groups",
        "Decarboxylation",
        "Hydrolysis",
        "Hydrogenation",
        "Oxidation of alcohol",
        "Nitrile to amide",
    ]

    def is_significant_reactant(smiles):
        """Determine if a reactant is significant (not a common reagent)"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        # Skip very small molecules (likely reagents)
        if mol.GetNumAtoms() < 3:
            return False

        # Skip common reagents
        for reagent in common_reagents:
            if reagent in smiles:
                return False

        return True

    def dfs_traverse(node, depth=0):
        nonlocal linear_steps_count, convergent_steps_count, total_steps

        if node["type"] == "reaction":
            total_steps += 1
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a known convergent reaction type
                is_convergent_rxn = False
                is_linear_rxn = False

                for rxn_type in convergent_reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_convergent_rxn = True
                        print(f"Found convergent reaction type: {rxn_type} at depth {depth}")
                        break

                if not is_convergent_rxn:
                    for rxn_type in linear_reaction_types:
                        if checker.check_reaction(rxn_type, rsmi):
                            is_linear_rxn = True
                            print(f"Found linear reaction type: {rxn_type} at depth {depth}")
                            break

                # Count significant reactants
                significant_reactants = [r for r in reactants if is_significant_reactant(r)]

                # Classify the step
                if is_convergent_rxn or len(significant_reactants) > 1:
                    print(
                        f"Found convergent step at depth {depth} with {len(significant_reactants)} significant reactants"
                    )
                    convergent_steps_count += 1
                elif is_linear_rxn or len(significant_reactants) <= 1:
                    print(f"Found linear step at depth {depth}")
                    linear_steps_count += 1

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    print(f"Total steps: {total_steps}")
    print(f"Linear steps: {linear_steps_count}")
    print(f"Convergent steps: {convergent_steps_count}")

    # If no steps were found, return False
    if total_steps == 0:
        return False

    # Consider it a linear strategy if:
    # 1. More than 70% of steps are linear, or
    # 2. There are at least 3 linear steps and linear steps outnumber convergent steps by at least 2:1
    linear_ratio = linear_steps_count / total_steps if total_steps > 0 else 0
    return (linear_ratio >= 0.7) or (
        linear_steps_count >= 3 and linear_steps_count >= convergent_steps_count * 2
    )

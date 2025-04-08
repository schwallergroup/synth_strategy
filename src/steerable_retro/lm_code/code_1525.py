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

root_data = "/home/andres/Documents/steerable_retro/data"

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
    This function detects if the synthetic route follows a linear strategy
    without convergent steps (each reaction has only one non-reagent reactant).
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            # Check if this is a known convergent reaction type
            convergent_reactions = [
                # Coupling reactions
                "Suzuki coupling with boronic acids",
                "Suzuki coupling with boronic acids OTf",
                "Suzuki coupling with boronic esters",
                "Suzuki coupling with boronic esters OTf",
                "Suzuki coupling with sulfonic esters",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                "Sonogashira acetylene_aryl halide",
                "Sonogashira alkyne_aryl halide",
                "Sonogashira acetylene_aryl OTf",
                "Sonogashira alkyne_aryl OTf",
                "Negishi coupling",
                "Stille reaction_vinyl",
                "Stille reaction_aryl",
                "Stille reaction_benzyl",
                "Stille reaction_allyl",
                "Stille reaction_vinyl OTf",
                "Stille reaction_aryl OTf",
                "Heck terminal vinyl",
                "Oxidative Heck reaction",
                "Heck reaction with vinyl ester and amine",
                "Ugi reaction",
                "Hiyama-Denmark Coupling",
                "Kumada cross-coupling",
                "Aryllithium cross-coupling",
                # Amide formations
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acyl chloride with secondary amine to amide",
                "Carboxylic acid with primary amine to amide",
                "Ester with primary amine to amide",
                "Ester with secondary amine to amide",
                # Other convergent reactions
                "Mitsunobu esterification",
                "Mitsunobu aryl ether",
                "Diels-Alder",
                "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                "Urea synthesis via isocyanate and primary amine",
                "Urea synthesis via isocyanate and secondary amine",
                "Aldol condensation",
                "Wittig reaction with triphenylphosphorane",
                "Reductive amination with aldehyde",
                "Reductive amination with ketone",
            ]

            for rxn_type in convergent_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    print(f"Found convergent reaction: {rxn_type}")
                    is_linear = False
                    return

            # Count non-reagent reactants (complex molecules)
            complex_reactants = 0
            for r in reactants:
                mol = Chem.MolFromSmiles(r)
                if mol:
                    num_atoms = mol.GetNumHeavyAtoms()
                    # Consider a molecule complex if it has more than 8 heavy atoms
                    # and is not a common reagent
                    if num_atoms > 8:
                        # Check if it's not a common reagent
                        is_reagent = False

                        # Check for common organometallic reagents
                        organometallic_fgs = [
                            "Magnesium halide",
                            "Alkyl lithium",
                            "Zinc halide",
                            "Boronic acid",
                            "Boronic ester",
                            "Tin",
                            "Silane",
                        ]
                        if any(checker.check_fg(fg, r) for fg in organometallic_fgs):
                            is_reagent = True

                        # Check for common protecting groups
                        protecting_fgs = [
                            "TMS ether protective group",
                            "Silyl protective group",
                            "Acetal/Ketal",
                            "Boc",
                        ]
                        if any(checker.check_fg(fg, r) for fg in protecting_fgs):
                            is_reagent = True

                        # Check for common reagent patterns
                        reagent_patterns = [
                            "B(OH)2",
                            "B(O",
                            "MgCl",
                            "MgBr",
                            "MgI",
                            "ZnCl",
                            "Li",
                            "SnBu3",
                            "Si(",
                            "SiMe3",
                            "TMS",
                            "P(Ph)3",
                            "PPh3",
                            "P(O",
                            "N(Bu)4",
                            "DMAP",
                            "EDC",
                            "DCC",
                            "HOBt",
                            "HATU",
                            "PyBOP",
                            "TBTU",
                            "DIPEA",
                            "NEt3",
                            "DBU",
                            "LDA",
                            "NaH",
                            "KH",
                            "NaOH",
                            "KOH",
                            "Boc2O",
                            "Fmoc",
                            "Cbz",
                            "TBDMS",
                            "TBDPS",
                            "THP",
                            "Ac2O",
                        ]
                        if any(pattern in r for pattern in reagent_patterns):
                            is_reagent = True

                        # Check for small reagents with functional groups
                        reagent_fgs = [
                            "Triflate",
                            "Mesylate",
                            "Tosylate",
                            "Azide",
                            "Isocyanate",
                            "Diazo",
                            "Acyl halide",
                            "Sulfonyl halide",
                            "Anhydride",
                        ]
                        if num_atoms <= 15 and any(checker.check_fg(fg, r) for fg in reagent_fgs):
                            is_reagent = True

                        # Check for common solvents and small molecules
                        if num_atoms <= 12 and any(
                            solv in r.lower()
                            for solv in [
                                "thf",
                                "dmf",
                                "dmso",
                                "meoh",
                                "etoh",
                                "acoh",
                                "h2o",
                                "dcm",
                                "chcl3",
                                "ether",
                                "toluene",
                                "acetone",
                                "acn",
                                "mecn",
                            ]
                        ):
                            is_reagent = True

                        if not is_reagent:
                            complex_reactants += 1
                            print(f"Found complex reactant: {r} with {num_atoms} atoms")

            if complex_reactants > 1:
                print(f"Found convergent step with {complex_reactants} complex reactants")
                is_linear = False

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return is_linear

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
    This function detects esterification reactions (carboxylic acid → ester)
    or the reverse (ester → carboxylic acid) in a synthetic route.
    """
    esterification_detected = False

    def dfs_traverse(node):
        nonlocal esterification_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction: {rsmi}")

            # Check for direct esterification reactions using the checker function
            esterification_reactions = [
                "Esterification of Carboxylic Acids",
                "Transesterification",
                "O-alkylation of carboxylic acids with diazo compounds",
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                "Schotten-Baumann to ester",
                "Oxidative esterification of primary alcohols",
            ]

            # Check for reverse esterification (hydrolysis) reactions
            hydrolysis_reactions = [
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                "Ester saponification (methyl deprotection)",
                "Ester saponification (alkyl deprotection)",
            ]

            # Check for forward esterification reactions
            for rxn_name in esterification_reactions:
                if checker.check_reaction(rxn_name, rsmi):
                    print(f"Esterification reaction detected: {rxn_name}")
                    esterification_detected = True
                    return

            # Check for reverse esterification reactions (hydrolysis)
            for rxn_name in hydrolysis_reactions:
                if checker.check_reaction(rxn_name, rsmi):
                    print(f"Ester hydrolysis reaction detected: {rxn_name}")
                    esterification_detected = True
                    return

            # As a fallback, check for the functional group transformation manually
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for forward esterification (carboxylic acid + alcohol → ester)
                carboxylic_acid_in_reactants = any(
                    checker.check_fg("Carboxylic acid", reactant) for reactant in reactants
                )
                alcohol_in_reactants = any(
                    checker.check_fg("Primary alcohol", reactant)
                    or checker.check_fg("Secondary alcohol", reactant)
                    or checker.check_fg("Tertiary alcohol", reactant)
                    or checker.check_fg("Phenol", reactant)
                    for reactant in reactants
                )
                ester_in_product = checker.check_fg("Ester", product)

                # Check for reverse esterification (ester → carboxylic acid)
                ester_in_reactants = any(
                    checker.check_fg("Ester", reactant) for reactant in reactants
                )
                carboxylic_acid_in_product = checker.check_fg("Carboxylic acid", product)

                print(
                    f"Manual check - Forward: Carboxylic acid in reactants: {carboxylic_acid_in_reactants}, "
                    f"Alcohol in reactants: {alcohol_in_reactants}, Ester in product: {ester_in_product}"
                )
                print(
                    f"Manual check - Reverse: Ester in reactants: {ester_in_reactants}, "
                    f"Carboxylic acid in product: {carboxylic_acid_in_product}"
                )

                # Forward esterification check
                if carboxylic_acid_in_reactants and alcohol_in_reactants and ester_in_product:
                    print(f"Forward esterification detected through functional group analysis")
                    esterification_detected = True
                    return

                # Reverse esterification check (hydrolysis)
                if ester_in_reactants and carboxylic_acid_in_product:
                    print(
                        f"Reverse esterification (hydrolysis) detected through functional group analysis"
                    )
                    esterification_detected = True
                    return

                # Check for acyl halide esterification
                acyl_halide_in_reactants = any(
                    checker.check_fg("Acyl halide", reactant) for reactant in reactants
                )
                if acyl_halide_in_reactants and alcohol_in_reactants and ester_in_product:
                    print(f"Acyl halide esterification detected through functional group analysis")
                    esterification_detected = True
                    return

            except Exception as e:
                print(f"Error in manual esterification check: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Final result: esterification_detected = {esterification_detected}")
    return esterification_detected

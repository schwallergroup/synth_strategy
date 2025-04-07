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
    This function detects if the synthesis route involves an amide coupling reaction.
    """
    amide_coupling = False

    def dfs_traverse(node):
        nonlocal amide_coupling

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amide coupling reactions directly using the checker function
                amide_coupling_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with ammonia to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with ammonia to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                    "Aminolysis of esters",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                ]

                for reaction_type in amide_coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Amide coupling detected: {reaction_type} in reaction: {rsmi}")
                        amide_coupling = True
                        return

                # If no specific reaction type matched, check for amide formation manually
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    has_amide = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )

                    if has_amide:
                        # Check if reactants contain amine and carboxylic acid/derivative
                        has_amine = False
                        has_carbonyl = False

                        for reactant in reactants:
                            if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                                "Secondary amine", reactant
                            ):
                                has_amine = True
                                print(f"Found amine in reactant: {reactant}")

                            if (
                                checker.check_fg("Carboxylic acid", reactant)
                                or checker.check_fg("Acyl halide", reactant)
                                or checker.check_fg("Ester", reactant)
                                or checker.check_fg("Anhydride", reactant)
                            ):
                                has_carbonyl = True
                                print(f"Found carbonyl in reactant: {reactant}")

                        if has_amine and has_carbonyl:
                            print(f"Amide coupling detected through manual check: {rsmi}")
                            amide_coupling = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return amide_coupling

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
    This function detects amide bond formation from carboxylic acid derivatives and amines.
    """
    amide_formation_detected = False

    def dfs_traverse(node):
        nonlocal amide_formation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a known amide formation reaction
            amide_formation_reactions = [
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
                "Carboxylic acid to amide conversion",
            ]

            for reaction_type in amide_formation_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Detected amide formation reaction: {reaction_type}")
                    amide_formation_detected = True
                    return

            # If no specific reaction type matched, check for functional groups
            # Check for amide in product
            if (
                not checker.check_fg("Primary amide", product)
                and not checker.check_fg("Secondary amide", product)
                and not checker.check_fg("Tertiary amide", product)
            ):
                return

            print(f"Found amide in product: {product}")

            # Check for carboxylic acid derivatives in reactants
            carboxylic_acid_derivatives = ["Carboxylic acid", "Ester", "Acyl halide", "Anhydride"]

            carboxylic_acid_present = False
            for reactant in reactants:
                for derivative in carboxylic_acid_derivatives:
                    if checker.check_fg(derivative, reactant):
                        carboxylic_acid_present = True
                        print(f"Found carboxylic acid derivative ({derivative}): {reactant}")
                        break

            # Check for amine in reactants
            amine_present = False
            for reactant in reactants:
                if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                    "Secondary amine", reactant
                ):
                    amine_present = True
                    print(f"Found amine: {reactant}")
                    break

            if carboxylic_acid_present and amine_present:
                # Additional verification: try to match with common amide formation patterns
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Check for C(=O)N pattern in product
                    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
                    if product_mol.HasSubstructMatch(amide_pattern):
                        amide_formation_detected = True
                        print("Amide formation detected based on functional groups!")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return amide_formation_detected

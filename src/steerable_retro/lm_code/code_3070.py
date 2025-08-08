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
    This function detects a strategy involving modification of a pyrazole core.
    """
    has_pyrazole_modification = False

    def dfs_traverse(node, depth=0):
        nonlocal has_pyrazole_modification

        if node["type"] == "reaction":
            try:
                # Extract product and reactants
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product has pyrazole
                if checker.check_ring("pyrazole", product_smiles):
                    print(f"Found product with pyrazole at depth {depth}: {product_smiles}")

                    # Case 1: Pyrazole modification - pyrazole in both reactant and product
                    reactant_has_pyrazole = False
                    for reactant in reactants_smiles:
                        if checker.check_ring("pyrazole", reactant):
                            print(f"Found reactant with pyrazole at depth {depth}: {reactant}")
                            reactant_has_pyrazole = True

                            # Check for any reaction that could modify pyrazoles
                            modification_reactions = [
                                "Friedel-Crafts acylation",
                                "Friedel-Crafts alkylation",
                                "N-alkylation of primary amines with alkyl halides",
                                "N-alkylation of secondary amines with alkyl halides",
                                "Suzuki coupling with boronic acids",
                                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                                "{pyrazole}",
                                "Esterification of Carboxylic Acids",
                                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                                "Acylation of primary amines",
                                "Acylation of secondary amines",
                                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                                "Oxidation of aldehydes to carboxylic acids",
                                "Oxidation of alcohol to carboxylic acid",
                            ]

                            # Check if any specific reaction type is detected
                            reaction_detected = False
                            for reaction_type in modification_reactions:
                                if checker.check_reaction(reaction_type, rsmi):
                                    print(
                                        f"Detected pyrazole modification reaction at depth {depth}: {reaction_type}"
                                    )
                                    reaction_detected = True
                                    has_pyrazole_modification = True
                                    break

                            # If no specific reaction detected but pyrazole is modified, consider it a valid modification
                            if not reaction_detected:
                                print(f"Detected generic pyrazole modification at depth {depth}")
                                has_pyrazole_modification = True

                    # Case 2: Pyrazole formation - pyrazole in product but not in reactants
                    if not reactant_has_pyrazole:
                        formation_reactions = [
                            "{pyrazole}",
                            "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                            "Huisgen 1,3 dipolar cycloaddition",
                            "Pyrazole formation",
                        ]

                        # Check if any specific formation reaction is detected
                        reaction_detected = False
                        for reaction_type in formation_reactions:
                            if checker.check_reaction(reaction_type, rsmi):
                                print(
                                    f"Detected pyrazole formation reaction at depth {depth}: {reaction_type}"
                                )
                                reaction_detected = True
                                has_pyrazole_modification = True
                                break

                        # If no specific reaction detected but pyrazole is formed, consider it a valid formation
                        if not reaction_detected:
                            print(f"Detected generic pyrazole formation at depth {depth}")
                            has_pyrazole_modification = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Process children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: {has_pyrazole_modification}")

    return has_pyrazole_modification

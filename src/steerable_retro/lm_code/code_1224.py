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
    Detects if the synthetic route involves a convergent synthesis with amide coupling.
    """
    amide_coupling_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an amide coupling reaction
                is_amide_coupling = False

                # Check for amide coupling reactions using the reaction checker
                amide_coupling_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Carboxylic acid with primary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                ]

                for reaction_name in amide_coupling_reactions:
                    if checker.check_reaction(reaction_name, rsmi):
                        is_amide_coupling = True
                        break

                # If not detected by reaction type, check for functional group transformation
                if not is_amide_coupling:
                    has_acid = False
                    has_amine = False
                    has_amide_product = False

                    for reactant in reactants:
                        try:
                            if checker.check_fg("Carboxylic acid", reactant):
                                has_acid = True
                            if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                                "Secondary amine", reactant
                            ):
                                has_amine = True
                        except:
                            continue

                    try:
                        if (
                            checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                        ):
                            has_amide_product = True
                    except:
                        pass

                    if has_acid and has_amine and has_amide_product and len(reactants) >= 2:
                        is_amide_coupling = True

                if is_amide_coupling:
                    # Check if this is a convergent step (multiple complex fragments coming together)
                    complex_reactants = 0
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if (
                                mol and mol.GetNumAtoms() > 6
                            ):  # Consider reactants with more than 6 atoms as complex
                                complex_reactants += 1
                        except:
                            continue

                    if complex_reactants >= 2:  # At least two complex fragments coming together
                        amide_coupling_depths.append(depth)
                        print(f"Found convergent amide coupling at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if there's at least one amide coupling at depth 1-5 (mid-synthesis)
    result = any(1 <= depth <= 5 for depth in amide_coupling_depths)
    print(f"Convergent synthesis with amide coupling detected: {result}")
    return result

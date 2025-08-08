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
    This function detects if the synthesis follows a linear strategy (no convergent steps).
    A linear synthesis builds the target molecule by adding one fragment at a time,
    while convergent synthesis combines multiple complex fragments.
    """
    is_linear = True

    # List of common solvents and reagents that shouldn't count as complex reactants
    common_reagents = ["Acetic acid", "Methanol", "Ethanol", "Water", "Triethylamine"]

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            parts = rsmi.split(">")
            if len(parts) >= 3:  # Ensure valid reaction SMILES
                reactants_smiles = parts[0].split(".")
                product_smiles = parts[2]

                # Skip simple reactions that are typically linear
                simple_reaction_types = [
                    "Acylation",
                    "Alkylation",
                    "Protection",
                    "Deprotection",
                    "Reduction",
                    "Oxidation",
                    "Hydrolysis",
                ]

                if any(
                    checker.check_reaction(rxn_type, rsmi) for rxn_type in simple_reaction_types
                ):
                    print(f"Simple reaction detected: {rsmi}")
                    # These reactions are typically linear even with multiple reactants
                    pass
                else:
                    # Count significant reactants (more than 7 atoms and not a common reagent)
                    significant_reactants = []

                    for r_smiles in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r_smiles)
                        if r_mol and r_mol.GetNumAtoms() > 7:
                            # Check if it's not a common reagent
                            is_reagent = False
                            for reagent in common_reagents:
                                if checker.check_fg(reagent, r_smiles):
                                    is_reagent = True
                                    break

                            if not is_reagent:
                                significant_reactants.append(r_smiles)

                    # If more than one significant reactant, analyze atom mapping to confirm convergence
                    if len(significant_reactants) > 1:
                        # Check if this is a coupling reaction (true convergent step)
                        is_coupling = any(
                            checker.check_reaction(rxn, rsmi)
                            for rxn in [
                                "Suzuki",
                                "Stille",
                                "Negishi",
                                "Heck",
                                "Sonogashira",
                                "Buchwald-Hartwig",
                                "Ullmann",
                            ]
                        )

                        if is_coupling or len(significant_reactants) > 2:
                            print(f"Convergent step detected in reaction: {rsmi}")
                            print(f"Significant reactants: {len(significant_reactants)}")
                            is_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return is_linear

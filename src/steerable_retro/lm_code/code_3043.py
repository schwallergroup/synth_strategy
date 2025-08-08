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
    This function detects a late-stage ketone reduction strategy where a ketone (C=O)
    is reduced to a methylene (CH2) group in one of the final synthetic steps.
    """
    found_ketone_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_ketone_reduction

        if node["type"] == "reaction" and depth <= 3:  # Focus on late-stage reactions
            print(f"Examining reaction at depth {depth}")
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Reaction SMILES: {rsmi}")
                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check for ketone in reactants
                ketone_found = False
                for reactant_smiles in reactants:
                    if checker.check_fg("Ketone", reactant_smiles):
                        ketone_found = True
                        ketone_reactant = reactant_smiles
                        print(f"Found ketone in reactant: {reactant_smiles}")
                        break

                if ketone_found:
                    # Check for various ketone reduction reactions
                    if checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi):
                        print("Detected 'Reduction of ketone to secondary alcohol' reaction")
                        if checker.check_fg("Secondary alcohol", product) and not checker.check_fg(
                            "Ketone", product
                        ):
                            print(f"Found late-stage ketone reduction to alcohol at depth {depth}")
                            found_ketone_reduction = True

                    elif checker.check_reaction(
                        "Reduction of aldehydes and ketones to alcohols", rsmi
                    ):
                        print("Detected 'Reduction of aldehydes and ketones to alcohols' reaction")
                        if checker.check_fg("Secondary alcohol", product) and not checker.check_fg(
                            "Ketone", product
                        ):
                            print(f"Found late-stage ketone reduction to alcohol at depth {depth}")
                            found_ketone_reduction = True

                    # Check for complete reduction reactions (ketone to methylene)
                    elif checker.check_reaction(
                        "Wolff-Kishner reduction", rsmi
                    ) or checker.check_reaction("Clemmensen reduction", rsmi):
                        print("Detected Wolff-Kishner or Clemmensen reduction")
                        if not checker.check_fg("Ketone", product) and not checker.check_fg(
                            "Secondary alcohol", product
                        ):
                            print(f"Found late-stage complete ketone reduction at depth {depth}")
                            found_ketone_reduction = True

                    # Generic reduction check
                    elif not checker.check_fg("Ketone", product):
                        # If ketone is gone but no specific reaction was identified
                        if checker.check_fg("Secondary alcohol", product):
                            print(
                                f"Found generic late-stage ketone reduction to alcohol at depth {depth}"
                            )
                            found_ketone_reduction = True
                        elif not checker.check_fg("Secondary alcohol", product):
                            # Look for other indicators of reduction
                            if (
                                "reduction" in rsmi.lower()
                                or "hydride" in rsmi.lower()
                                or "H2" in rsmi
                            ):
                                print(
                                    f"Found generic late-stage complete ketone reduction at depth {depth}"
                                )
                                found_ketone_reduction = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {found_ketone_reduction}")
    return found_ketone_reduction

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
    This function detects if the synthesis involves a late-stage esterification
    with a complex fragment containing a trifluoromethoxy group.
    """
    esterification_detected = False

    def is_complex_fragment(smiles):
        """Check if a molecule is a complex fragment (has multiple rings or functional groups)"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Check for complexity: rings
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() >= 1:
            print(f"Complex due to rings: {ring_info.NumRings()} rings")
            return True

        # Count different functional groups
        fg_count = 0
        for fg in [
            "Aromatic halide",
            "Ether",
            "Ester",
            "Amide",
            "Nitrile",
            "Nitro group",
            "Carboxylic acid",
            "Phenol",
            "Alcohol",
            "Trifluoro group",
        ]:
            if checker.check_fg(fg, smiles):
                fg_count += 1
                print(f"Found functional group: {fg}")

        if fg_count >= 2:  # At least 2 functional groups for complexity
            print(f"Complex due to functional groups: {fg_count}")
            return True

        # If molecule has at least 12 atoms, consider it complex
        if mol.GetNumHeavyAtoms() >= 12:
            print(f"Complex due to size: {mol.GetNumHeavyAtoms()} heavy atoms")
            return True

        return False

    def has_trifluoromethoxy(smiles):
        """Check if a molecule contains a trifluoromethoxy group"""
        # Use checker to detect trifluoro group
        if checker.check_fg("Trifluoro group", smiles):
            print(f"Found trifluoro group in: {smiles}")

            # Check if it's connected to oxygen (trifluoromethoxy)
            if "OC(F)(F)F" in smiles or "OCF3" in smiles or "OC(F)F" in smiles:
                print(f"Found trifluoromethoxy pattern in SMILES: {smiles}")
                return True

            # Check using ether functional group
            if checker.check_fg("Ether", smiles):
                print(f"Found ether with trifluoro group: {smiles}")
                return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal esterification_detected

        if node["type"] == "reaction" and depth <= 2:  # Late stage (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check for esterification reactions - try multiple types
                is_esterification = False
                for rxn_type in [
                    "Esterification of Carboxylic Acids",
                    "Transesterification",
                    "Schotten-Baumann to ester",
                    "O-alkylation of carboxylic acids with diazo compounds",
                    "Oxidative esterification of primary alcohols",
                ]:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Found esterification reaction type: {rxn_type}")
                        is_esterification = True
                        break

                # Fallback check: look for ester formation
                if not is_esterification:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]
                    if any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    ) and checker.check_fg("Ester", product):
                        print("Found esterification based on functional group transformation")
                        is_esterification = True

                if is_esterification:
                    # Extract reactants
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if any reactant has trifluoromethoxy group and is complex
                    for r in reactants:
                        if has_trifluoromethoxy(r):
                            print(f"Found reactant with trifluoromethoxy group: {r}")
                            if is_complex_fragment(r):
                                print(
                                    f"Confirmed complex fragment with trifluoromethoxy group: {r}"
                                )
                                esterification_detected = True
                                return  # Early return once found

                    # Also check if the product contains a trifluoromethoxy group
                    # This handles cases where the trifluoromethoxy group is in the product
                    if has_trifluoromethoxy(product):
                        print(f"Found product with trifluoromethoxy group: {product}")
                        # Check if any reactant is complex
                        for r in reactants:
                            if is_complex_fragment(r):
                                print(f"Confirmed complex fragment in reactants: {r}")
                                esterification_detected = True
                                return  # Early return once found

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Late-stage esterification with complex fragment: {esterification_detected}")
    return esterification_detected

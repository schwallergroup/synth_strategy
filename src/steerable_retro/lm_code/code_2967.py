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
    This function detects if the synthetic route involves late-stage protection (in the final 3 steps).
    """
    late_protection_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal late_protection_detected

        if node["type"] == "reaction" and depth <= 2:  # Final 3 steps (depth 0, 1, 2)
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    print(f"Analyzing reaction at depth {depth}: {rsmi}")

                    # Check if this is a protection reaction
                    protection_reactions = [
                        "Alcohol protection with silyl ethers",
                        "Boc amine protection",
                        "Boc amine protection explicit",
                        "Boc amine protection with Boc anhydride",
                        "Boc amine protection (ethyl Boc)",
                        "Boc amine protection of secondary amine",
                        "Boc amine protection of primary amine",
                        "Protection of carboxylic acid",
                        "Schotten-Baumann to ester",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    ]

                    for reaction_type in protection_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(
                                f"Late-stage protection detected: {reaction_type} at depth {depth}"
                            )
                            late_protection_detected = True
                            return

                    # Check for protection groups that might not be covered by reaction types
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if not product_mol:
                        return

                    # List of protection groups to check
                    protection_groups = [
                        "TMS ether protective group",
                        "Silyl protective group",
                        "Acetal/Ketal",
                        "Boc",
                        "Carbamic ester",  # Cbz, Fmoc, etc.
                    ]

                    # Check if protection group is formed in this reaction
                    for pg in protection_groups:
                        # Check if protection group is in product
                        if checker.check_fg(pg, product_smiles):
                            # Check if protection group is absent from ALL reactants
                            if not any(
                                checker.check_fg(pg, reactant) for reactant in reactants_smiles
                            ):
                                # Verify this is likely a protection reaction by checking for relevant functional groups
                                if (
                                    pg == "TMS ether protective group"
                                    or pg == "Silyl protective group"
                                ):
                                    for reactant in reactants_smiles:
                                        if (
                                            checker.check_fg("Primary alcohol", reactant)
                                            or checker.check_fg("Secondary alcohol", reactant)
                                            or checker.check_fg("Tertiary alcohol", reactant)
                                        ):
                                            print(
                                                f"Late-stage protection detected: {pg} added to alcohol at depth {depth}"
                                            )
                                            late_protection_detected = True
                                            return
                                elif pg == "Acetal/Ketal":
                                    for reactant in reactants_smiles:
                                        if checker.check_fg(
                                            "Aldehyde", reactant
                                        ) or checker.check_fg("Ketone", reactant):
                                            print(
                                                f"Late-stage protection detected: {pg} formed from carbonyl at depth {depth}"
                                            )
                                            late_protection_detected = True
                                            return
                                elif pg == "Boc":
                                    for reactant in reactants_smiles:
                                        if checker.check_fg(
                                            "Primary amine", reactant
                                        ) or checker.check_fg("Secondary amine", reactant):
                                            print(
                                                f"Late-stage protection detected: {pg} added to amine at depth {depth}"
                                            )
                                            late_protection_detected = True
                                            return
                                elif pg == "Carbamic ester":  # Cbz, Fmoc protection
                                    for reactant in reactants_smiles:
                                        if checker.check_fg(
                                            "Primary amine", reactant
                                        ) or checker.check_fg("Secondary amine", reactant):
                                            print(
                                                f"Late-stage protection detected: Cbz/carbamate protection of amine at depth {depth}"
                                            )
                                            late_protection_detected = True
                                            return

                    # Special check for Cbz protection (benzyloxycarbonyl)
                    if "OCc1ccccc1" in rsmi and any(
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                        for reactant in reactants_smiles
                    ):
                        # Check if a carbamate is formed
                        if "O=C(OCc1ccccc1)N" in rsmi or "O=C(OCc1ccccc1)[N" in rsmi:
                            print(
                                f"Late-stage protection detected: Cbz protection at depth {depth}"
                            )
                            late_protection_detected = True
                            return
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {late_protection_detected}")
    return late_protection_detected

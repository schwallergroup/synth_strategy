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
    Detects if the synthesis involves amide bond formation using a protection/deprotection strategy.
    """
    # Track reactions and their depths
    amide_formations = []  # List of (depth, rsmi) tuples
    protections = []  # List of (depth, rsmi, protected_group) tuples
    deprotections = []  # List of (depth, rsmi, deprotected_group) tuples

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amide formation reactions
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with primary amine to amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Schotten-Baumann to ester",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    "Acyl chloride with ammonia to amide",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with ammonia to amide",
                    "Ester with secondary amine to amide",
                    "{Schotten-Baumann_amide}",
                ]

                for reaction_type in amide_formation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected amide formation at depth {depth}: {reaction_type}")
                        amide_formations.append((depth, rsmi))
                        break

                # If no specific reaction type matched, check for amide formation by structure
                if not any(af[1] == rsmi for af in amide_formations):
                    # Check for amine in reactants
                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                    )
                    # Check for carboxylic acid or derivative in reactants
                    has_acid_derivative = any(
                        checker.check_fg("Carboxylic acid", r)
                        or checker.check_fg("Acyl halide", r)
                        or checker.check_fg("Ester", r)
                        or checker.check_fg("Anhydride", r)
                        for r in reactants
                    )
                    # Check for amide in product
                    has_amide = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )

                    if has_amine and has_acid_derivative and has_amide:
                        print(f"Detected amide formation by structure analysis at depth {depth}")
                        amide_formations.append((depth, rsmi))

                # Check for protection reactions
                protection_reactions = {
                    "Protection of carboxylic acid": "Carboxylic acid",
                    "Alcohol protection with silyl ethers": "Primary alcohol",
                    "Boc amine protection": "Primary amine",
                    "Boc amine protection explicit": "Primary amine",
                    "Boc amine protection with Boc anhydride": "Primary amine",
                    "Boc amine protection (ethyl Boc)": "Primary amine",
                    "Boc amine protection of secondary amine": "Secondary amine",
                    "Boc amine protection of primary amine": "Primary amine",
                    "TMS ether protective group": "Primary alcohol",
                    "Silyl protective group": "Primary alcohol",
                    "Acetal/Ketal": "Aldehyde",
                    "Boc": "Primary amine",
                }

                for reaction_type, protected_group in protection_reactions.items():
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Detected protection at depth {depth}: {reaction_type} - {protected_group}"
                        )
                        protections.append((depth, rsmi, protected_group))
                        break

                # Additional check for protection by functional group analysis
                if not any(p[1] == rsmi for p in protections):
                    # Check for Boc group in product
                    if any(checker.check_fg("Boc", r) for r in reactants) and checker.check_fg(
                        "Boc", product
                    ):
                        if any(checker.check_fg("Primary amine", r) for r in reactants) and not any(
                            checker.check_fg("Primary amine", product)
                        ):
                            print(f"Detected Boc protection by structure analysis at depth {depth}")
                            protections.append((depth, rsmi, "Primary amine"))

                    # Check for silyl protection
                    if any(
                        checker.check_fg("Silyl protective group", r) for r in reactants
                    ) and checker.check_fg("Silyl protective group", product):
                        if any(
                            checker.check_fg("Primary alcohol", r) for r in reactants
                        ) and not any(checker.check_fg("Primary alcohol", product)):
                            print(
                                f"Detected silyl protection by structure analysis at depth {depth}"
                            )
                            protections.append((depth, rsmi, "Primary alcohol"))

                # Check for deprotection reactions
                deprotection_reactions = {
                    "Deprotection of carboxylic acid": "Carboxylic acid",
                    "Alcohol deprotection from silyl ethers": "Primary alcohol",
                    "Alcohol deprotection from silyl ethers (double)": "Primary alcohol",
                    "Alcohol deprotection from silyl ethers (diol)": "Primary alcohol",
                    "Boc amine deprotection": "Primary amine",
                    "Boc amine deprotection of guanidine": "Primary amine",
                    "Boc amine deprotection to NH-NH2": "Primary amine",
                    "Tert-butyl deprotection of amine": "Primary amine",
                    "Hydroxyl benzyl deprotection": "Primary alcohol",
                    "Carboxyl benzyl deprotection": "Carboxylic acid",
                    "COOH ethyl deprotection": "Carboxylic acid",
                    "N-glutarimide deprotection": "Primary amine",
                    "Phthalimide deprotection": "Primary amine",
                    "TMS deprotection from alkyne": "Alkyne",
                    "Acetal hydrolysis to aldehyde": "Aldehyde",
                    "Ketal hydrolysis to ketone": "Ketone",
                    "Ester saponification (methyl deprotection)": "Carboxylic acid",
                    "Ester saponification (alkyl deprotection)": "Carboxylic acid",
                    "Cleavage of methoxy ethers to alcohols": "Primary alcohol",
                    "Cleavage of alkoxy ethers to alcohols": "Primary alcohol",
                    "Ether cleavage to primary alcohol": "Primary alcohol",
                }

                for reaction_type, deprotected_group in deprotection_reactions.items():
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Detected deprotection at depth {depth}: {reaction_type} - {deprotected_group}"
                        )
                        deprotections.append((depth, rsmi, deprotected_group))
                        break

                # Additional check for deprotection by functional group analysis
                if not any(d[1] == rsmi for d in deprotections):
                    # Check for Boc group removal
                    if any(checker.check_fg("Boc", r) for r in reactants) and not checker.check_fg(
                        "Boc", product
                    ):
                        if not any(
                            checker.check_fg("Primary amine", r) for r in reactants
                        ) and checker.check_fg("Primary amine", product):
                            print(
                                f"Detected Boc deprotection by structure analysis at depth {depth}"
                            )
                            deprotections.append((depth, rsmi, "Primary amine"))

                    # Check for silyl deprotection
                    if any(
                        checker.check_fg("Silyl protective group", r) for r in reactants
                    ) and not checker.check_fg("Silyl protective group", product):
                        if not any(
                            checker.check_fg("Primary alcohol", r) for r in reactants
                        ) and checker.check_fg("Primary alcohol", product):
                            print(
                                f"Detected silyl deprotection by structure analysis at depth {depth}"
                            )
                            deprotections.append((depth, rsmi, "Primary alcohol"))

                    # Check for ester hydrolysis to carboxylic acid (deprotection)
                    if any(
                        checker.check_fg("Ester", r) for r in reactants
                    ) and not checker.check_fg("Ester", product):
                        if checker.check_fg("Carboxylic acid", product):
                            print(
                                f"Detected ester hydrolysis (deprotection) by structure analysis at depth {depth}"
                            )
                            deprotections.append((depth, rsmi, "Carboxylic acid"))

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have amide formation
    if not amide_formations:
        print("No amide formation detected")
        return False

    # Check if we have protection or deprotection
    if not protections and not deprotections:
        print("No protection or deprotection detected")
        return False

    # Get the depths
    amide_depths = [d for d, _ in amide_formations]
    protection_depths = [d for d, _, _ in protections]
    deprotection_depths = [d for d, _, _ in deprotections]

    print(f"Amide formation depths: {amide_depths}")
    print(f"Protection depths: {protection_depths}")
    print(f"Deprotection depths: {deprotection_depths}")

    # Check if there's a protection before or at the same depth as amide formation (higher or equal depth in retrosynthesis)
    protection_before_amide = any(pd >= ad for pd in protection_depths for ad in amide_depths)

    # Check if there's a deprotection after or at the same depth as amide formation (lower or equal depth in retrosynthesis)
    deprotection_after_amide = any(dd <= ad for dd in deprotection_depths for ad in amide_depths)

    # Check if the protected/deprotected groups are relevant to amide formation
    relevant_groups = [
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Carboxylic acid",
        "Acyl halide",
        "Ester",
        "Anhydride",
        "Primary amide",
        "Secondary amide",
        "Tertiary amide",
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Aldehyde",
        "Ketone",
    ]

    protected_groups = [pg for _, _, pg in protections]
    deprotected_groups = [dg for _, _, dg in deprotections]

    print(f"Protected groups: {protected_groups}")
    print(f"Deprotected groups: {deprotected_groups}")

    relevant_protection = any(pg in relevant_groups for pg in protected_groups)
    relevant_deprotection = any(dg in relevant_groups for dg in deprotected_groups)

    # Final result: We need amide formation AND either:
    # 1. Protection before/at amide formation with relevant group, OR
    # 2. Deprotection after/at amide formation with relevant group
    result = len(amide_formations) > 0 and (
        (protection_before_amide and relevant_protection)
        or (deprotection_after_amide and relevant_deprotection)
    )

    if result:
        print("Detected amide formation with protection/deprotection strategy")
    else:
        if amide_formations:
            print("Amide formation detected but no relevant protection/deprotection strategy found")
            print(
                f"Protection before amide: {protection_before_amide}, Relevant protection: {relevant_protection}"
            )
            print(
                f"Deprotection after amide: {deprotection_after_amide}, Relevant deprotection: {relevant_deprotection}"
            )
        elif protections or deprotections:
            print("Protection/deprotection detected but no amide formation found")
        else:
            print("Neither amide formation nor protection/deprotection detected")

    return result

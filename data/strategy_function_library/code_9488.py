from typing import Tuple, Dict, List
import copy
from rdkit.Chem import AllChem, rdFMCS
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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

from pathlib import Path
root_data = Path(__file__).parent.parent

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


PROTECTING_GROUPS_OF_INTEREST = [
    "Sulfonamide",
    "Boc",
    "TMS ether protective group",
    "Silyl protective group",
    "Acetal/Ketal",
]

def main(route) -> Tuple[bool, Dict]:
    """ 
    This function detects a synthetic strategy involving the exchange of specific protecting groups, including Sulfonamide, Boc, TMS ether protective group, Silyl protective group, and Acetal/Ketal. It identifies both direct, single-step exchanges (PG1 -> PG2) and sequential, multi-step exchanges where a deprotection of one group is followed by a protection with another group in a subsequent step.
    """
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Initialize tracking variables
    protecting_group_reactions = {}  # Dictionary to store reactions by protecting group type
    direct_exchanges = []  # Store direct exchanges between protecting groups

    # Define protecting groups to track
    protecting_groups = PROTECTING_GROUPS_OF_INTEREST

    # Initialize dictionaries
    for pg in protecting_groups:
        protecting_group_reactions[pg] = []  # Store (depth, is_protection) tuples

    def dfs_traverse(node, depth=0):
        nonlocal findings_json, direct_exchanges, protecting_group_reactions
        if node["type"] == "mol":
            # This logic is too broad and a source of false positives, so it is not used.
            pass

        elif node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants = parts[0].split(".")
            products = parts[2].split(".")

            if not reactants or not products:
                return

            # Check for free amines
            reactant_has_free_amine = False
            for r in reactants:
                if r:
                    if checker.check_fg("Primary amine", r):
                        reactant_has_free_amine = True
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    if checker.check_fg("Secondary amine", r):
                        reactant_has_free_amine = True
                        if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

            product_has_free_amine = False
            for p in products:
                if p:
                    if checker.check_fg("Primary amine", p):
                        product_has_free_amine = True
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    if checker.check_fg("Secondary amine", p):
                        product_has_free_amine = True
                        if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

            # Check for direct exchanges between protecting groups
            for pg1 in protecting_groups:
                reactant_has_pg1 = any(checker.check_fg(pg1, r) for r in reactants if r)
                if reactant_has_pg1 and pg1 not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append(pg1)

                product_has_pg1 = any(checker.check_fg(pg1, p) for p in products if p)
                if product_has_pg1 and pg1 not in findings_json["atomic_checks"]["functional_groups"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append(pg1)

                for pg2 in protecting_groups:
                    if pg1 == pg2:
                        continue

                    reactant_has_pg2 = any(checker.check_fg(pg2, r) for r in reactants if r)
                    if reactant_has_pg2 and pg2 not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(pg2)

                    product_has_pg2 = any(checker.check_fg(pg2, p) for p in products if p)
                    if product_has_pg2 and pg2 not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(pg2)

                    # Direct exchange: pg1 in reactants, pg2 in products, no pg1 in products, no pg2 in reactants
                    if (
                        reactant_has_pg1
                        and product_has_pg2
                        and not product_has_pg1
                        and not reactant_has_pg2
                    ):
                        direct_exchanges.append((depth, pg1, pg2))
                        print(f"Found direct {pg1}-to-{pg2} exchange at depth {depth}")
                        print(f"Reaction SMILES: {rsmi}")
                        # Add structural constraint for direct exchange
                        constraint = {
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "protecting_group_deprotection",
                                    "protecting_group_protection"
                                ],
                                "condition": "Events must occur within the same reaction step and involve two different protecting groups from the specified list."
                            }
                        }
                        if constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(constraint)

            # Check for Sulfonamide reactions
            sulfonamide_deprotection_found = checker.check_reaction("Sulfonamide deprotection", rsmi)
            if sulfonamide_deprotection_found:
                if "Sulfonamide deprotection" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Sulfonamide deprotection")
            
            reactant_has_sulfonamide = any(checker.check_fg("Sulfonamide", r) for r in reactants if r)
            if reactant_has_sulfonamide and "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")

            product_has_sulfonamide = any(checker.check_fg("Sulfonamide", p) for p in products if p)
            if product_has_sulfonamide and "Sulfonamide" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Sulfonamide")

            if sulfonamide_deprotection_found or (
                reactant_has_sulfonamide
                and product_has_free_amine
                and not product_has_sulfonamide
            ):
                protecting_group_reactions["Sulfonamide"].append(
                    (depth, False)
                )  # False = deprotection
                print(f"Found sulfonamide deprotection at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")

            sulfonamide_prot_schotten_baumann_primary = checker.check_reaction("Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi)
            if sulfonamide_prot_schotten_baumann_primary:
                if "Sulfonamide synthesis (Schotten-Baumann) primary amine" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Sulfonamide synthesis (Schotten-Baumann) primary amine")

            sulfonamide_prot_schotten_baumann_secondary = checker.check_reaction("Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi)
            if sulfonamide_prot_schotten_baumann_secondary:
                if "Sulfonamide synthesis (Schotten-Baumann) secondary amine" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Sulfonamide synthesis (Schotten-Baumann) secondary amine")

            if (
                sulfonamide_prot_schotten_baumann_primary
                or sulfonamide_prot_schotten_baumann_secondary
            ) or (
                reactant_has_free_amine
                and product_has_sulfonamide
                and not reactant_has_sulfonamide
            ):
                protecting_group_reactions["Sulfonamide"].append((depth, True))  # True = protection
                print(f"Found sulfonamide protection at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")

            # Check for Boc reactions
            boc_amine_protection = checker.check_reaction("Boc amine protection", rsmi)
            if boc_amine_protection:
                if "Boc amine protection" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Boc amine protection")
            boc_amine_protection_explicit = checker.check_reaction("Boc amine protection explicit", rsmi)
            if boc_amine_protection_explicit:
                if "Boc amine protection explicit" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Boc amine protection explicit")
            boc_amine_protection_anhydride = checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
            if boc_amine_protection_anhydride:
                if "Boc amine protection with Boc anhydride" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Boc amine protection with Boc anhydride")
            boc_amine_protection_ethyl = checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
            if boc_amine_protection_ethyl:
                if "Boc amine protection (ethyl Boc)" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Boc amine protection (ethyl Boc)")
            boc_amine_protection_secondary = checker.check_reaction("Boc amine protection of secondary amine", rsmi)
            if boc_amine_protection_secondary:
                if "Boc amine protection of secondary amine" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Boc amine protection of secondary amine")
            boc_amine_protection_primary = checker.check_reaction("Boc amine protection of primary amine", rsmi)
            if boc_amine_protection_primary:
                if "Boc amine protection of primary amine" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Boc amine protection of primary amine")

            reactant_has_boc = any(checker.check_fg("Boc", r) for r in reactants if r)
            if reactant_has_boc and "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Boc")
            product_has_boc = any(checker.check_fg("Boc", p) for p in products if p)
            if product_has_boc and "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Boc")

            if (
                boc_amine_protection
                or boc_amine_protection_explicit
                or boc_amine_protection_anhydride
                or boc_amine_protection_ethyl
                or boc_amine_protection_secondary
                or boc_amine_protection_primary
            ) or (
                reactant_has_free_amine
                and product_has_boc
                and not reactant_has_boc
            ):
                protecting_group_reactions["Boc"].append((depth, True))  # True = protection
                print(f"Found Boc protection at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")

            boc_amine_deprotection = checker.check_reaction("Boc amine deprotection", rsmi)
            if boc_amine_deprotection:
                if "Boc amine deprotection" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Boc amine deprotection")
            boc_amine_deprotection_guanidine = checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
            if boc_amine_deprotection_guanidine:
                if "Boc amine deprotection of guanidine" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Boc amine deprotection of guanidine")
            boc_amine_deprotection_nh2 = checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
            if boc_amine_deprotection_nh2:
                if "Boc amine deprotection to NH-NH2" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Boc amine deprotection to NH-NH2")
            tert_butyl_deprotection = checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
            if tert_butyl_deprotection:
                if "Tert-butyl deprotection of amine" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Tert-butyl deprotection of amine")

            if (
                boc_amine_deprotection
                or boc_amine_deprotection_guanidine
                or boc_amine_deprotection_nh2
                or tert_butyl_deprotection
            ) or (
                reactant_has_boc
                and product_has_free_amine
                and not product_has_boc
            ):
                protecting_group_reactions["Boc"].append((depth, False))  # False = deprotection
                print(f"Found Boc deprotection at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")

            # Check for silyl protection/deprotection
            silyl_protection = checker.check_reaction("Alcohol protection with silyl ethers", rsmi)
            if silyl_protection:
                if "Alcohol protection with silyl ethers" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Alcohol protection with silyl ethers")
            if silyl_protection:
                protecting_group_reactions["Silyl protective group"].append((depth, True))
                print(f"Found silyl protection at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")
                if "Silyl protective group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Silyl protective group")

            silyl_deprotection = checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
            if silyl_deprotection:
                if "Alcohol deprotection from silyl ethers" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Alcohol deprotection from silyl ethers")
            silyl_deprotection_double = checker.check_reaction("Alcohol deprotection from silyl ethers (double)", rsmi)
            if silyl_deprotection_double:
                if "Alcohol deprotection from silyl ethers (double)" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Alcohol deprotection from silyl ethers (double)")
            silyl_deprotection_diol = checker.check_reaction("Alcohol deprotection from silyl ethers (diol)", rsmi)
            if silyl_deprotection_diol:
                if "Alcohol deprotection from silyl ethers (diol)" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Alcohol deprotection from silyl ethers (diol)")

            if (
                silyl_deprotection
                or silyl_deprotection_double
                or silyl_deprotection_diol
            ):
                protecting_group_reactions["Silyl protective group"].append((depth, False))
                print(f"Found silyl deprotection at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")
                if "Silyl protective group" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Silyl protective group")

            # Check for acetal/ketal protection/deprotection
            acetal_ketal_protection_aldehyde_ketone = checker.check_reaction("Aldehyde or ketone acetalization", rsmi)
            if acetal_ketal_protection_aldehyde_ketone:
                if "Aldehyde or ketone acetalization" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Aldehyde or ketone acetalization")
            acetal_ketal_protection_diol = checker.check_reaction("Diol acetalization", rsmi)
            if acetal_ketal_protection_diol:
                if "Diol acetalization" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Diol acetalization")

            if acetal_ketal_protection_aldehyde_ketone or acetal_ketal_protection_diol:
                protecting_group_reactions["Acetal/Ketal"].append((depth, True))
                print(f"Found acetal/ketal protection at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")
                if "Acetal/Ketal" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Acetal/Ketal")

            acetal_hydrolysis_diol = checker.check_reaction("Acetal hydrolysis to diol", rsmi)
            if acetal_hydrolysis_diol:
                if "Acetal hydrolysis to diol" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Acetal hydrolysis to diol")
            acetal_hydrolysis_aldehyde = checker.check_reaction("Acetal hydrolysis to aldehyde", rsmi)
            if acetal_hydrolysis_aldehyde:
                if "Acetal hydrolysis to aldehyde" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Acetal hydrolysis to aldehyde")
            ketal_hydrolysis_ketone = checker.check_reaction("Ketal hydrolysis to ketone", rsmi)
            if ketal_hydrolysis_ketone:
                if "Ketal hydrolysis to ketone" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Ketal hydrolysis to ketone")

            if (
                acetal_hydrolysis_diol
                or acetal_hydrolysis_aldehyde
                or ketal_hydrolysis_ketone
            ):
                protecting_group_reactions["Acetal/Ketal"].append((depth, False))
                print(f"Found acetal/ketal deprotection at depth {depth}")
                print(f"Reaction SMILES: {rsmi}")
                if "Acetal/Ketal" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Acetal/Ketal")

        # Process children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "mol"
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Analyze results to determine if protecting group exchange strategy is present
    has_protecting_group_exchange = False

    # Check for direct exchange
    if direct_exchanges:
        print(f"Direct protecting group exchanges found: {direct_exchanges}")
        has_protecting_group_exchange = True

    # Check for sequential deprotection/protection patterns
    for pg1 in protecting_groups:
        deprotections_pg1 = [
            (depth, is_prot) for depth, is_prot in protecting_group_reactions[pg1] if not is_prot
        ]

        for pg2 in protecting_groups:
            if pg1 == pg2:
                continue

            protections_pg2 = [
                (depth, is_prot) for depth, is_prot in protecting_group_reactions[pg2] if is_prot
            ]

            print(f"{pg1} deprotections: {deprotections_pg1}")
            print(f"{pg2} protections: {protections_pg2}")

            # In a forward synthesis, deprotection of PG1 occurs at an earlier step (higher retro depth)
            # than the protection with PG2 (lower retro depth).
            # Therefore, we must check if p_depth < d_depth.
            for d_depth, _ in deprotections_pg1:
                for p_depth, _ in protections_pg2:
                    if p_depth < d_depth:
                        print(
                            f"Found protecting group exchange pattern: {pg1} deprotection at depth {d_depth} followed by {pg2} protection at depth {p_depth}"
                        )
                        has_protecting_group_exchange = True
                        # Add structural constraint for sequential exchange
                        constraint = {
                            "type": "sequence",
                            "details": {
                                "ordered_events": [
                                    "protecting_group_deprotection",
                                    "protecting_group_protection"
                                ],
                                "condition": "The deprotection of one protecting group must occur in an earlier synthetic step than the protection with a different protecting group."
                            }
                        }
                        if constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(constraint)

    print(f"Protecting group exchange detected: {has_protecting_group_exchange}")
    return has_protecting_group_exchange, findings_json

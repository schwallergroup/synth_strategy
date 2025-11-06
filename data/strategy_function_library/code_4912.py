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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

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


PROTECTION_PAIRS = [
    # (protection_reaction, deprotection_reaction, protected_fg, original_fg)
    ("Boc amine protection", "Boc amine deprotection", "Carbamic ester", "Primary amine"),
    ("Boc amine protection", "Boc amine deprotection", "Carbamic ester", "Secondary amine"),
    ("Alcohol protection with silyl ethers", "Alcohol deprotection from silyl ethers", "Silyl protective group", "Primary alcohol"),
    ("Alcohol protection with silyl ethers", "Alcohol deprotection from silyl ethers", "Silyl protective group", "Secondary alcohol"),
    ("Alcohol protection with silyl ethers", "Alcohol deprotection from silyl ethers", "Silyl protective group", "Tertiary alcohol"),
    ("Protection of carboxylic acid", "Deprotection of carboxylic acid", "Ester", "Carboxylic acid"),
    ("Aldehyde or ketone acetalization", "Acetal/Ketal", "Acetal/Ketal", "Aldehyde"),
    ("Aldehyde or ketone acetalization", "Ketal hydrolysis to ketone", "Acetal/Ketal", "Ketone"),
    ("Diol acetalization", "Acetal hydrolysis to diol", "Acetal/Ketal", "Primary alcohol"),
]

# Define the full strategy JSON for reference
STRATEGY_JSON = {
  "function_id": "code_4912",
  "filepath": "../data/merged_good_perf/code_4912.py",
  "description": "This function detects a protection-deprotection sequence strategy where functional groups are protected and later deprotected in the synthesis.",
  "atomic_checks": {
    "named_reactions": [
      "Acetal hydrolysis to aldehyde",
      "Acetal hydrolysis to diol",
      "Alcohol deprotection from silyl ethers",
      "Alcohol protection with silyl ethers",
      "Aldehyde or ketone acetalization",
      "Boc amine deprotection",
      "Boc amine protection",
      "Deprotection of carboxylic acid",
      "Diol acetalization",
      "Ketal hydrolysis to ketone",
      "Protection of carboxylic acid"
    ],
    "ring_systems": [],
    "functional_groups": [
      "Acetal/Ketal",
      "Aldehyde",
      "Carbamic ester",
      "Carboxylic acid",
      "Ester",
      "Ketone",
      "Primary alcohol",
      "Primary amine",
      "Secondary alcohol",
      "Secondary amine",
      "Silyl protective group",
      "Tertiary alcohol"
    ]
  },
  "structural_constraints": [
    {
      "type": "sequence",
      "details": {
        "before": "Boc amine protection",
        "after": "Boc amine deprotection"
      }
    },
    {
      "type": "sequence",
      "details": {
        "before": "Alcohol protection with silyl ethers",
        "after": "Alcohol deprotection from silyl ethers"
      }
    },
    {
      "type": "sequence",
      "details": {
        "before": "Protection of carboxylic acid",
        "after": "Deprotection of carboxylic acid"
      }
    },
    {
      "type": "sequence",
      "details": {
        "before": "Aldehyde or ketone acetalization",
        "after": "Acetal hydrolysis to aldehyde"
      }
    },
    {
      "type": "sequence",
      "details": {
        "before": "Aldehyde or ketone acetalization",
        "after": "Ketal hydrolysis to ketone"
      }
    },
    {
      "type": "sequence",
      "details": {
        "before": "Diol acetalization",
        "after": "Acetal hydrolysis to diol"
      }
    }
  ]
}

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a protection-deprotection sequence strategy where
    functional groups are protected and later deprotected in the synthesis.
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
    result = False

    # Track protection and deprotection events with molecule identifiers
    protection_events = []
    deprotection_events = []

    def dfs_traverse(node, depth=0, molecule_path=None):
        nonlocal findings_json
        if molecule_path is None:
            molecule_path = []

        if node["type"] == "reaction":
            try:
                # Extract reactants and product from reaction SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check for protection reactions
                for prot_rxn, deprot_rxn, protected_fg, original_fg in PROTECTION_PAIRS:
                    # Check for protection event
                    if checker.check_reaction(prot_rxn, rsmi):
                        if prot_rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(prot_rxn)
                        # Verify the functional group transformation
                        if checker.check_fg(original_fg, reactants_part):
                            if original_fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(original_fg)
                        if checker.check_fg(protected_fg, product_part):
                            if protected_fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(protected_fg)

                        if checker.check_fg(original_fg, reactants_part) and checker.check_fg(
                            protected_fg, product_part
                        ):
                            print(f"Found {prot_rxn} at depth {depth}: {rsmi}")
                            protection_events.append((prot_rxn, deprot_rxn, depth, rsmi))

                    # Check for deprotection event
                    if checker.check_reaction(deprot_rxn, rsmi):
                        if deprot_rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(deprot_rxn)
                        # Verify the functional group transformation
                        if checker.check_fg(protected_fg, reactants_part):
                            if protected_fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(protected_fg)
                        if checker.check_fg(original_fg, product_part):
                            if original_fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(original_fg)

                        if checker.check_fg(protected_fg, reactants_part) and checker.check_fg(
                            original_fg, product_part
                        ):
                            print(f"Found {deprot_rxn} at depth {depth}: {rsmi}")
                            deprotection_events.append((prot_rxn, deprot_rxn, depth, rsmi))

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a chemical node
            next_depth = depth + 1

        # Process children (moving deeper in the retrosynthetic tree)
        for child in node.get("children", []):
            dfs_traverse(child, next_depth, molecule_path + [node])

    # Start traversal
    dfs_traverse(route)

    # Check if we have both protection and deprotection events
    if protection_events and deprotection_events:
        print(f"Protection events: {protection_events}")
        print(f"Deprotection events: {deprotection_events}")

        # Check if the same type of protection/deprotection occurs in the correct sequence
        for prot_rxn, deprot_rxn, prot_depth, prot_rsmi in protection_events:
            for prot_rxn_check, deprot_rxn_check, deprot_depth, deprot_rsmi in deprotection_events:
                # Check if it's the same protection/deprotection pair
                if prot_rxn == prot_rxn_check and deprot_rxn == deprot_rxn_check:
                    # In retrosynthesis, protection is at higher depth (earlier in synthesis)
                    # than deprotection (later in synthesis)
                    print(
                        f"Comparing: Protection at depth {prot_depth}, deprotection at depth {deprot_depth}"
                    )
                    if prot_depth > deprot_depth:
                        print(f"Found protection-deprotection sequence: {prot_rxn} -> {deprot_rxn}")
                        print(
                            f"Protection at depth {prot_depth}, deprotection at depth {deprot_depth}"
                        )
                        result = True
                        # Add the corresponding structural constraint
                        for constraint in STRATEGY_JSON["structural_constraints"]:
                            if constraint["type"] == "sequence" and \
                               constraint["details"]["before"] == prot_rxn and \
                               constraint["details"]["after"] == deprot_rxn:
                                if constraint not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append(constraint)
                        # No need to continue checking other pairs if one is found
                        # as the function returns True for the first found sequence
                        return result, findings_json

    print("No protection-deprotection sequence found")
    return result, findings_json

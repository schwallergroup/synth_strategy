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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a protection-deprotection sequence strategy,
    specifically looking for Boc protection followed by deprotection,
    or other common protection-deprotection sequences.
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

    # Track protection and deprotection events with atom mapping
    protection_events = {
        "boc": [],  # [(depth, atom_maps)]
        "phthalimide": [],  # [(depth, atom_maps)]
        "silyl": [],  # [(depth, atom_maps)]
        "benzyl": [],  # [(depth, atom_maps)]
        "acetal": [],  # [(depth, atom_maps)]
    }

    deprotection_events = {
        "boc": [],  # [(depth, atom_maps)]
        "phthalimide": [],  # [(depth, atom_maps)]
        "silyl": [],  # [(depth, atom_maps)]
        "benzyl": [],  # [(depth, atom_maps)]
        "acetal": [],  # [(depth, atom_maps)]
    }

    def extract_atom_maps(mol_smiles):
        """Extract atom mapping numbers from a molecule SMILES"""
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return []

        atom_maps = []
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() > 0:
                atom_maps.append(atom.GetAtomMapNum())
        return atom_maps

    def dfs_traverse(node, depth=0):
        nonlocal protection_events, deprotection_events, findings_json
        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                products_part = rsmi.split(">")[-1]

                reactants = reactants_part.split(".")
                product = products_part

                # Extract atom maps for tracking
                reactant_maps = [extract_atom_maps(r) for r in reactants]
                product_maps = extract_atom_maps(product)

                # Check for Boc protection reactions
                boc_prot_reactions = [
                    "Boc amine protection",
                    "Boc amine protection explicit",
                    "Boc amine protection with Boc anhydride",
                    "Boc amine protection (ethyl Boc)",
                    "Boc amine protection of secondary amine",
                    "Boc amine protection of primary amine"
                ]
                for r_name in boc_prot_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        print(f"Found {r_name} at depth {depth}")
                        protection_events["boc"].append((depth, product_maps))
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                # Check for Boc deprotection reactions
                boc_deprot_reactions = [
                    "Boc amine deprotection",
                    "Boc amine deprotection of guanidine",
                    "Boc amine deprotection to NH-NH2",
                    "Tert-butyl deprotection of amine"
                ]
                for r_name in boc_deprot_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        print(f"Found {r_name} at depth {depth}")
                        deprotection_events["boc"].append((depth, product_maps))
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                # Check for Phthalimide protection/deprotection
                phthalimide_prot_reactions = [
                    "Phthalic anhydride to phthalimide",
                    "Phthalimide_formation"
                ]
                phthalimide_prot_found = False
                for r_name in phthalimide_prot_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        print(f"Found {r_name} at depth {depth}")
                        protection_events["phthalimide"].append((depth, product_maps))
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        phthalimide_prot_found = True
                        break
                
                if not phthalimide_prot_found and (checker.check_fg("Phthalimide", product) and not any(checker.check_fg("Phthalimide", r) for r in reactants)):
                    print(f"Found Phthalimide protection (FG check) at depth {depth}")
                    protection_events["phthalimide"].append((depth, product_maps))
                    if "Phthalimide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Phthalimide")

                phthalimide_deprot_reaction = "Phthalimide deprotection"
                if checker.check_reaction(phthalimide_deprot_reaction, rsmi):
                    print(f"Found {phthalimide_deprot_reaction} at depth {depth}")
                    deprotection_events["phthalimide"].append((depth, product_maps))
                    if phthalimide_deprot_reaction not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(phthalimide_deprot_reaction)

                # Check for silyl protection/deprotection
                silyl_prot_reaction = "Alcohol protection with silyl ethers"
                if checker.check_reaction(silyl_prot_reaction, rsmi):
                    print(f"Found {silyl_prot_reaction} at depth {depth}")
                    protection_events["silyl"].append((depth, product_maps))
                    if silyl_prot_reaction not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(silyl_prot_reaction)

                silyl_deprot_reactions = [
                    "Alcohol deprotection from silyl ethers",
                    "Alcohol deprotection from silyl ethers (double)",
                    "Alcohol deprotection from silyl ethers (diol)",
                    "TMS deprotection from alkyne"
                ]
                for r_name in silyl_deprot_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        print(f"Found {r_name} at depth {depth}")
                        deprotection_events["silyl"].append((depth, product_maps))
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                # Check for benzyl deprotection
                benzyl_deprot_reactions = [
                    "Hydroxyl benzyl deprotection",
                    "Carboxyl benzyl deprotection"
                ]
                for r_name in benzyl_deprot_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        print(f"Found {r_name} at depth {depth}")
                        deprotection_events["benzyl"].append((depth, product_maps))
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                # Check for acetal protection/deprotection
                acetal_prot_reactions = [
                    "Aldehyde or ketone acetalization",
                    "Diol acetalization"
                ]
                for r_name in acetal_prot_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        print(f"Found {r_name} at depth {depth}")
                        protection_events["acetal"].append((depth, product_maps))
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

                acetal_deprot_reactions = [
                    "Acetal hydrolysis to diol",
                    "Acetal hydrolysis to aldehyde",
                    "Ketal hydrolysis to ketone"
                ]
                for r_name in acetal_deprot_reactions:
                    if checker.check_reaction(r_name, rsmi):
                        print(f"Found {r_name} at depth {depth}")
                        deprotection_events["acetal"].append((depth, product_maps))
                        if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(r_name)

        # Process children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            # Depth remains the same when traversing from reaction to chemical node
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' type or any other type that should increase depth
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have valid protection-deprotection sequences
    sequence_found = False

    for protect_type in protection_events:
        # Case 1: Both protection and deprotection events found
        if protection_events[protect_type] and deprotection_events[protect_type]:
            print(f"Found {protect_type} protection and deprotection events")

            # Check if there's at least one protection event at a higher depth than a deprotection event
            for prot_depth, prot_maps in protection_events[protect_type]:
                for deprot_depth, deprot_maps in deprotection_events[protect_type]:
                    # In retrosynthesis, higher depth = earlier in synthesis
                    # We want protection to be earlier in synthesis than deprotection
                    if prot_depth > deprot_depth:
                        # Check if there's overlap in atom maps to ensure it's the same functional group
                        # If atom mapping is incomplete, we'll assume it's a valid sequence
                        if (
                            not prot_maps
                            or not deprot_maps
                            or any(m in deprot_maps for m in prot_maps)
                        ):
                            print(
                                f"Valid {protect_type} protection-deprotection sequence found: protection at depth {prot_depth}, deprotection at depth {deprot_depth}"
                            )
                            sequence_found = True
                            # Add structural constraint for sequence
                            if {"type": "sequence", "details": {"description": "A protection reaction for a specific chemical group (e.g., Boc, Silyl, Acetal) must occur at a greater depth (earlier in the synthesis) than a deprotection reaction for the same group.", "antecedent_event_group": ["Boc amine protection", "Boc amine protection explicit", "Boc amine protection with Boc anhydride", "Boc amine protection (ethyl Boc)", "Boc amine protection of secondary amine", "Boc amine protection of primary amine", "Phthalic anhydride to phthalimide", "Phthalimide_formation", "Alcohol protection with silyl ethers", "Aldehyde or ketone acetalization", "Diol acetalization"], "consequent_event_group": ["Boc amine deprotection", "Boc amine deprotection of guanidine", "Boc amine deprotection to NH-NH2", "Tert-butyl deprotection of amine", "Phthalimide deprotection", "Alcohol deprotection from silyl ethers", "Alcohol deprotection from silyl ethers (double)", "Alcohol deprotection from silyl ethers (diol)", "TMS deprotection from alkyne", "Hydroxyl benzyl deprotection", "Carboxyl benzyl deprotection", "Acetal hydrolysis to diol", "Acetal hydrolysis to aldehyde", "Ketal hydrolysis to ketone"], "condition": "depth(antecedent) > depth(consequent) for events of the same protecting group family"}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "sequence", "details": {"description": "A protection reaction for a specific chemical group (e.g., Boc, Silyl, Acetal) must occur at a greater depth (earlier in the synthesis) than a deprotection reaction for the same group.", "antecedent_event_group": ["Boc amine protection", "Boc amine protection explicit", "Boc amine protection with Boc anhydride", "Boc amine protection (ethyl Boc)", "Boc amine protection of secondary amine", "Boc amine protection of primary amine", "Phthalic anhydride to phthalimide", "Phthalimide_formation", "Alcohol protection with silyl ethers", "Aldehyde or ketone acetalization", "Diol acetalization"], "consequent_event_group": ["Boc amine deprotection", "Boc amine deprotection of guanidine", "Boc amine deprotection to NH-NH2", "Tert-butyl deprotection of amine", "Phthalimide deprotection", "Alcohol deprotection from silyl ethers", "Alcohol deprotection from silyl ethers (double)", "Alcohol deprotection from silyl ethers (diol)", "TMS deprotection from alkyne", "Hydroxyl benzyl deprotection", "Carboxyl benzyl deprotection", "Acetal hydrolysis to diol", "Acetal hydrolysis to aldehyde", "Ketal hydrolysis to ketone"], "condition": "depth(antecedent) > depth(consequent) for events of the same protecting group family"}})

        # Case 2: Only deprotection found - this is still valid as protection might be implicit
        elif deprotection_events[protect_type]:
            print(f"Found {protect_type} deprotection events without explicit protection")
            # Common protecting groups that might be implicitly added
            common_protecting_groups = [
                "boc",
                "phthalimide",
                "silyl",
                "benzyl",
                "acetal",
            ]
            if protect_type in common_protecting_groups:
                print(f"Assuming implicit {protect_type} protection - valid sequence")
                sequence_found = True
                # Add structural constraint for co-occurrence
                if {"type": "co-occurrence", "details": {"description": "Alternatively, the strategy is considered valid if at least one deprotection reaction for a common protecting group is found, assuming the protection step was implicit in the starting material.", "targets": ["Boc amine deprotection", "Boc amine deprotection of guanidine", "Boc amine deprotection to NH-NH2", "Tert-butyl deprotection of amine", "Phthalimide deprotection", "Alcohol deprotection from silyl ethers", "Alcohol deprotection from silyl ethers (double)", "Alcohol deprotection from silyl ethers (diol)", "TMS deprotection from alkyne", "Hydroxyl benzyl deprotection", "Carboxyl benzyl deprotection", "Acetal hydrolysis to diol", "Acetal hydrolysis to aldehyde", "Ketal hydrolysis to ketone"], "min_occurrences": 1, "logic": "OR"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"description": "Alternatively, the strategy is considered valid if at least one deprotection reaction for a common protecting group is found, assuming the protection step was implicit in the starting material.", "targets": ["Boc amine deprotection", "Boc amine deprotection of guanidine", "Boc amine deprotection to NH-NH2", "Tert-butyl deprotection of amine", "Phthalimide deprotection", "Alcohol deprotection from silyl ethers", "Alcohol deprotection from silyl ethers (double)", "Alcohol deprotection from silyl ethers (diol)", "TMS deprotection from alkyne", "Hydroxyl benzyl deprotection", "Carboxyl benzyl deprotection", "Acetal hydrolysis to diol", "Acetal hydrolysis to aldehyde", "Ketal hydrolysis to ketone"], "min_occurrences": 1, "logic": "OR"}})

    # Print summary
    for protect_type in protection_events:
        print(
            f"{protect_type.capitalize()} protection found: {len(protection_events[protect_type]) > 0}, {protect_type} deprotection found: {len(deprotection_events[protect_type]) > 0}"
        )

    return sequence_found, findings_json

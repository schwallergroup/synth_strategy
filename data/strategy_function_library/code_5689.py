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


AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with secondary amine to amide",
    "Ester with ammonia to amide",
    "Ester with secondary amine to amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
]

AROMATIC_SCAFFOLDS_OF_INTEREST = [
    "benzene",
    "pyridine",
    "pyrrole",
    "furan",
    "thiophene",
    "imidazole",
    "pyrazole",
    "oxazole",
    "thiazole",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
]

AMIDE_TYPES = ["Primary amide", "Secondary amide", "Tertiary amide"]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy involving a sequence of nitrogen oxidation state changes (nitro → amine → nitrile → amine → amide) while maintaining a core aromatic scaffold. The specific named reactions for amide formation and the aromatic scaffolds being searched for are defined in the `AMIDE_FORMATION_REACTIONS` and `AROMATIC_SCAFFOLDS_OF_INTEREST` lists, respectively.
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

    # Track nitrogen oxidation state changes
    n_oxidation_changes = 0

    # Track the sequence of transformations
    transformation_sequence = []

    # Track functional groups present at each step
    nitro_present = False
    amine_present = False
    nitrile_present = False
    amide_present = False
    aromatic_present = False

    def dfs_traverse(node, depth):
        nonlocal n_oxidation_changes, nitro_present, amine_present, nitrile_present, amide_present, aromatic_present, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")
            product = product_part

            # Check for nitro reduction to amine
            if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                print(f"Detected nitro reduction at depth {depth}")
                n_oxidation_changes += 1
                transformation_sequence.append("nitro_to_amine")
                findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
                findings_json["atomic_checks"]["named_reactions"].append("nitro_to_amine")
            # Fallback check for nitro reduction
            elif any(checker.check_fg("Nitro group", r) for r in reactants) and checker.check_fg(
                "Primary amine", product
            ):
                if not any(checker.check_fg("Primary amine", r) for r in reactants):
                    print(f"Detected nitro to amine conversion (fallback) at depth {depth}")
                    n_oxidation_changes += 1
                    transformation_sequence.append("nitro_to_amine")
                    findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
                    findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    findings_json["atomic_checks"]["named_reactions"].append("nitro_to_amine")

            # Check for amine to nitrile conversion
            if any(checker.check_fg("Primary amine", r) for r in reactants) and checker.check_fg(
                "Nitrile", product
            ):
                if not any(checker.check_fg("Nitrile", r) for r in reactants):
                    print(f"Detected amine to nitrile conversion at depth {depth}")
                    n_oxidation_changes += 1
                    transformation_sequence.append("amine_to_nitrile")
                    findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                    findings_json["atomic_checks"]["named_reactions"].append("amine_to_nitrile")

            # Check for nitrile to amine reduction
            if checker.check_reaction("Reduction of nitrile to amine", rsmi):
                print(f"Detected nitrile to amine reduction (reaction) at depth {depth}")
                n_oxidation_changes += 1
                transformation_sequence.append("nitrile_to_amine")
                findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitrile to amine")
                findings_json["atomic_checks"]["named_reactions"].append("nitrile_to_amine")
            elif any(checker.check_fg("Nitrile", r) for r in reactants) and checker.check_fg(
                "Primary amine", product
            ):
                if not any(checker.check_fg("Primary amine", r) for r in reactants):
                    print(f"Detected nitrile to amine reduction at depth {depth}")
                    n_oxidation_changes += 1
                    transformation_sequence.append("nitrile_to_amine")
                    findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                    findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    findings_json["atomic_checks"]["named_reactions"].append("nitrile_to_amine")

            # Check for amine to amide conversion
            amide_rxn_found = False
            for rxn in AMIDE_FORMATION_REACTIONS:
                if checker.check_reaction(rxn, rsmi):
                    amide_rxn_found = True
                    findings_json["atomic_checks"]["named_reactions"].append(rxn)
            if amide_rxn_found:
                print(f"Detected amine to amide conversion at depth {depth}")
                n_oxidation_changes += 1
                transformation_sequence.append("amine_to_amide")
                findings_json["atomic_checks"]["named_reactions"].append("amine_to_amide")
            # Fallback check for amide formation
            elif any(checker.check_fg("Primary amine", r) for r in reactants) and any(
                checker.check_fg(amide_type, product) for amide_type in AMIDE_TYPES
            ):
                if not any(
                    any(checker.check_fg(amide_type, r) for amide_type in AMIDE_TYPES)
                    for r in reactants
                ):
                    print(f"Detected amine to amide conversion (fallback) at depth {depth}")
                    n_oxidation_changes += 1
                    transformation_sequence.append("amine_to_amide")
                    findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    for amide_type in AMIDE_TYPES:
                        if checker.check_fg(amide_type, product):
                            findings_json["atomic_checks"]["functional_groups"].append(amide_type)
                    findings_json["atomic_checks"]["named_reactions"].append("amine_to_amide")

        elif node["type"] == "mol":
            # Check for functional groups in molecules
            mol_smiles = node["smiles"]

            if checker.check_fg("Nitro group", mol_smiles):
                nitro_present = True
                findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
            if checker.check_fg("Primary amine", mol_smiles):
                amine_present = True
                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
            if checker.check_fg("Nitrile", mol_smiles):
                nitrile_present = True
                findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
            for amide_type in AMIDE_TYPES:
                if checker.check_fg(amide_type, mol_smiles):
                    amide_present = True
                    findings_json["atomic_checks"]["functional_groups"].append(amide_type)

            # Check for aromatic scaffold
            for ring in AROMATIC_SCAFFOLDS_OF_INTEREST:
                if checker.check_ring(ring, mol_smiles):
                    aromatic_present = True
                    findings_json["atomic_checks"]["ring_systems"].append(ring)

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'mol' or other types
                dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route, 0)

    # Check if we have the required sequence in the correct order
    print(f"Transformation sequence: {transformation_sequence}")

    # In retrosynthetic traversal, the sequence appears in reverse order
    # nitro → amine → nitrile → amine → amide becomes amide → amine → nitrile → amine → nitro
    sequence_pattern = ["amine_to_amide", "nitrile_to_amine", "amine_to_nitrile", "nitro_to_amine"]

    # Check if the transformations appear in the correct relative order
    sequence_present = False

    # Check if all required transformations are present
    all_transformations_present = all(t in transformation_sequence for t in sequence_pattern)

    # Check if they appear in the correct relative order
    if all_transformations_present:
        # Get the indices of each transformation in the sequence
        indices = {t: transformation_sequence.index(t) for t in sequence_pattern}
        # Check if the order is correct (in retrosynthetic direction)
        if (
            indices["amine_to_amide"] < indices["nitrile_to_amine"]
            and indices["nitrile_to_amine"] < indices["nitro_to_amine"]
        ):
            # The amine_to_nitrile should be somewhere between nitrile_to_amine and nitro_to_amine
            if (
                indices["amine_to_nitrile"] > indices["nitrile_to_amine"]
                and indices["amine_to_nitrile"] < indices["nitro_to_amine"]
            ):
                sequence_present = True

    # Alternative check for sequence - more flexible
    # If the strict check failed, try a more relaxed check
    if not sequence_present:
        # Check if we have at least 3 of the 4 key transformations
        key_transformations = [
            "nitro_to_amine",
            "amine_to_nitrile",
            "nitrile_to_amine",
            "amine_to_amide",
        ]
        present_key_transformations = [
            t for t in key_transformations if t in transformation_sequence
        ]

        if len(present_key_transformations) >= 3:
            # Check if they appear in roughly the correct order
            if (
                "nitro_to_amine" in transformation_sequence
                and "nitrile_to_amine" in transformation_sequence
                and "amine_to_amide" in transformation_sequence
            ):
                nitro_idx = transformation_sequence.index("nitro_to_amine")
                nitrile_idx = transformation_sequence.index("nitrile_to_amine")
                amide_idx = transformation_sequence.index("amine_to_amide")

                # In retrosynthetic direction: amide → amine → nitro
                if nitro_idx > nitrile_idx > amide_idx:
                    print("Detected correct general sequence order (relaxed check)")
                    sequence_present = True
                    findings_json["structural_constraints"].append({
                        "type": "sequence",
                        "details": {
                            "ordered_events": [
                                "amine_to_amide",
                                "nitrile_to_amine",
                                "nitro_to_amine"
                            ],
                            "description": "The route must contain a specific sequence of transformations. The check is relaxed, requiring at least these three events to occur in the specified retrosynthetic order (e.g., amide formation occurs in an earlier step than nitrile reduction, which occurs earlier than nitro reduction)."
                        }
                    })

    print(f"Total nitrogen oxidation state changes: {n_oxidation_changes}")
    print(
        f"Functional groups present: Nitro: {nitro_present}, Amine: {amine_present}, "
        f"Nitrile: {nitrile_present}, Amide: {amide_present}, Aromatic: {aromatic_present}"
    )
    print(f"Sequence present: {sequence_present}")

    result = n_oxidation_changes >= 3 and sequence_present and aromatic_present

    if n_oxidation_changes >= 3:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "nitrogen_oxidation_state_change",
                "operator": ">=",
                "value": 3,
                "description": "The route must contain at least 3 key transformations involving nitrogen oxidation state changes (nitro_to_amine, amine_to_nitrile, nitrile_to_amine, amine_to_amide)."
            }
        })
    if aromatic_present:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "any_AROMATIC_SCAFFOLDS_OF_INTEREST"
                ],
                "description": "At least one of the specified aromatic scaffolds must be present as a molecule anywhere in the synthesis route."
            }
        })

    # Ensure unique entries in atomic_checks lists
    for key in findings_json["atomic_checks"]:
        findings_json["atomic_checks"][key] = list(set(findings_json["atomic_checks"][key]))

    return result, findings_json

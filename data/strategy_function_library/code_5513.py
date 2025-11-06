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


# Refactored lists for enumeration
N_HETEROCYCLES = [
    "pyrazole", "imidazole", "triazole", "tetrazole", "pyrrole", "indole",
    "pyridine", "piperidine", "piperazine", "morpholine", "pyrrolidine",
    "azetidine", "aziridine", "azepane", "diazepane", "benzimidazole",
    "benzotriazole", "indazole", "purine",
]

ALKYLATING_LEAVING_GROUPS = [
    "Primary halide", "Secondary halide", "Tertiary halide",
    "Mesylate", "Tosylate", "Triflate",
]

# Corrected and refactored list of N-alkylation reactions
N_ALKYLATION_REACTIONS = [
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides",
    "Methylation",
    "N-methylation",
    "Eschweiler-Clarke Primary Amine Methylation",
    "Eschweiler-Clarke Secondary Amine Methylation",
    "Reductive methylation of primary amine with formaldehyde",
    "Alkylation of amines",
    "Mitsunobu reaction",
    "Mitsunobu_imide",
    "Mitsunobu_phenole",
    "Mitsunobu_sulfonamide",
    "Mitsunobu_tetrazole_1",
    "Mitsunobu_tetrazole_2",
    "Mitsunobu_tetrazole_3",
    "Mitsunobu_tetrazole_4",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage N-alkylation with a complex, amine-containing fragment.

    This strategy is identified in the final two steps of a synthesis if:
    1. The reaction is a known N-alkylation type (from `N_ALKYLATION_REACTIONS`).
    2. One reactant is a nitrogen heterocycle with an alkylatable N-H (from `N_HETEROCYCLES`).
    3. Another reactant contains a leaving group (from `ALKYLATING_LEAVING_GROUPS`), an amine, and has at least 5 atoms.
    4. The product confirms the transformation by retaining the heterocycle and amine moieties while losing the leaving group.
    """
    late_stage_alkylation_detected = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_alkylation_detected, findings_json

        if node["type"] == "reaction" and depth <= 2:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                is_n_alkylation = False
                for rxn in N_ALKYLATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        is_n_alkylation = True
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                
                if is_n_alkylation:
                    # Structural Constraint: positional - within_last_three_stages
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "N-alkylation reaction",
                            "position": "within_last_three_stages"
                        }
                    })

                    has_n_heterocycle = False
                    heterocycle_reactant = None
                    heterocycle_type = None

                    for reactant in reactants_smiles:
                        for heterocycle in N_HETEROCYCLES:
                            if checker.check_ring(heterocycle, reactant):
                                if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                                if (
                                    checker.check_fg("Secondary amine", reactant)
                                    or checker.check_fg("Primary amine", reactant)
                                    or "[nH]" in reactant
                                ):
                                    if checker.check_fg("Secondary amine", reactant) and "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                                    if checker.check_fg("Primary amine", reactant) and "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                        findings_json["atomic_checks"]["functional_groups"].append("Primary amine")

                                    has_n_heterocycle = True
                                    heterocycle_reactant = reactant
                                    heterocycle_type = heterocycle
                                    break
                        if has_n_heterocycle:
                            break

                    alkylating_agent = None
                    leaving_group_type = None

                    for reactant in reactants_smiles:
                        if reactant == heterocycle_reactant:
                            continue

                        for lg in ALKYLATING_LEAVING_GROUPS:
                            if checker.check_fg(lg, reactant):
                                if lg not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append(lg)
                                alkylating_agent = reactant
                                leaving_group_type = lg
                                break
                        if alkylating_agent:
                            break

                    has_complex_amine_fragment = False
                    if alkylating_agent:
                        has_amine = (
                            checker.check_fg("Primary amine", alkylating_agent)
                            or checker.check_fg("Secondary amine", alkylating_agent)
                            or checker.check_fg("Tertiary amine", alkylating_agent)
                            or checker.check_fg("Aniline", alkylating_agent)
                        )
                        if checker.check_fg("Primary amine", alkylating_agent) and "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        if checker.check_fg("Secondary amine", alkylating_agent) and "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                        if checker.check_fg("Tertiary amine", alkylating_agent) and "Tertiary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")
                        if checker.check_fg("Aniline", alkylating_agent) and "Aniline" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aniline")

                        mol = Chem.MolFromSmiles(alkylating_agent)
                        atom_count = mol.GetNumAtoms() if mol else 0

                        if has_amine and atom_count >= 5:
                            has_complex_amine_fragment = True
                            # Structural Constraint: count - atoms_in_alkylating_agent
                            findings_json["structural_constraints"].append({
                                "type": "count",
                                "details": {
                                    "target": "atoms_in_alkylating_agent",
                                    "operator": ">=",
                                    "value": 5,
                                    "description": "The reactant serving as the alkylating agent must contain at least 5 atoms."
                                }
                            })

                    if has_n_heterocycle and has_complex_amine_fragment:
                        # Structural Constraint: co-occurrence - reactants
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "scope": "reactants",
                                "targets": [
                                    "N-heterocycle with alkylatable nitrogen",
                                    "alkylating agent with amine group"
                                ],
                                "description": "The reaction must have one reactant that is a specified N-heterocycle with a primary/secondary amine or [nH] group, and another reactant that contains a specified leaving group and an amine."
                            }
                        })

                        product_has_heterocycle = checker.check_ring(
                            heterocycle_type, product_smiles
                        )
                        if product_has_heterocycle and heterocycle_type not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(heterocycle_type)

                        product_has_amine = (
                            checker.check_fg("Primary amine", product_smiles)
                            or checker.check_fg("Secondary amine", product_smiles)
                            or checker.check_fg("Tertiary amine", product_smiles)
                            or checker.check_fg("Aniline", product_smiles)
                        )
                        if checker.check_fg("Primary amine", product_smiles) and "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                        if checker.check_fg("Secondary amine", product_smiles) and "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")
                        if checker.check_fg("Tertiary amine", product_smiles) and "Tertiary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary amine")
                        if checker.check_fg("Aniline", product_smiles) and "Aniline" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aniline")

                        if product_has_heterocycle and product_has_amine:
                            # Structural Constraint: co-occurrence - product
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "scope": "product",
                                    "targets": [
                                        "N-heterocycle",
                                        "amine_group"
                                    ],
                                    "description": "The product must contain both the original N-heterocycle ring system and an amine group."
                                }
                            })

                            product_has_leaving_group = False
                            if leaving_group_type:
                                product_has_leaving_group = checker.check_fg(
                                    leaving_group_type, product_smiles
                                )

                            if not product_has_leaving_group:
                                # Structural Constraint: negation - leaving_group
                                findings_json["structural_constraints"].append({
                                    "type": "negation",
                                    "details": {
                                        "scope": "product",
                                        "target": "leaving_group",
                                        "description": "The product must not contain the leaving group that was present on the alkylating agent reactant."
                                    }
                                })
                                late_stage_alkylation_detected = True
            except Exception:
                pass

        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_alkylation_detected, findings_json

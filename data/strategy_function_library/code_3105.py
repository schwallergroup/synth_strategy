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


HETEROCYCLES_OF_INTEREST = [
    "furan", "pyrrole", "thiophene", "pyridine", "pyrazole", "imidazole",
    "oxazole", "thiazole", "pyrimidine", "pyrazine", "pyridazine", "triazole",
    "tetrazole", "indole", "benzimidazole", "benzoxazole", "benzothiazole",
]

ACYLATION_REACTIONS = [
    "Friedel-Crafts acylation",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of secondary amines",
    "Acylation of primary amines",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
]

HALOGENATION_REACTIONS = [
    "Aromatic fluorination", "Aromatic chlorination", "Aromatic bromination",
    "Aromatic iodination", "Chlorination", "Fluorination", "Iodination", "Bromination",
]

ACYL_FGS = [
    "Ketone", "Primary amide", "Secondary amide", "Tertiary amide",
    "Ester", "Carboxylic acid", "Acyl halide",
]

HALOGEN_FGS = [
    "Aromatic halide", "Primary halide", "Secondary halide", "Tertiary halide",
    "Alkenyl halide", "Haloalkyne", "Triflate",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a multi-step strategy involving specific heterocycles. The strategy is flagged if either: (1) A heterocycle from HETEROCYCLES_OF_INTEREST is formed early in the synthesis, followed by acylation and then halogenation of the scaffold in later stages. (2) A starting material already contains a heterocycle from the list, and is subsequently acylated and then halogenated. Acylation and halogenation are identified by specific reaction names (see ACYLATION_REACTIONS and HALOGENATION_REACTIONS) or by the net addition of specific functional groups (see ACYL_FGS and HALOGEN_FGS).
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

    found_heterocycle = False
    found_acylation = False
    found_halogenation = False
    heterocycle_depth = -1
    acylation_depth = -1
    halogenation_depth = -1

    final_product_smiles = route["smiles"] if route["type"] == "mol" else ""

    final_product_heterocycle = None
    for heterocycle in HETEROCYCLES_OF_INTEREST:
        if checker.check_ring(heterocycle, final_product_smiles):
            final_product_heterocycle = heterocycle
            found_heterocycle = True
            findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
            break

    if not found_heterocycle:
        return False, findings_json

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycle, found_acylation, found_halogenation
        nonlocal heterocycle_depth, acylation_depth, halogenation_depth, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                if not found_heterocycle or heterocycle_depth == -1:
                    for heterocycle in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(heterocycle, product_smiles):
                            if not any(checker.check_ring(heterocycle, r) for r in reactants_smiles):
                                found_heterocycle = True
                                heterocycle_depth = depth
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                                break

                if not found_acylation:
                    for acylation_reaction in ACYLATION_REACTIONS:
                        if checker.check_reaction(acylation_reaction, rsmi):
                            if checker.check_ring(final_product_heterocycle, product_smiles):
                                found_acylation = True
                                acylation_depth = depth
                                findings_json["atomic_checks"]["named_reactions"].append(acylation_reaction)
                                break

                if not found_acylation:
                    for acyl_fg in ACYL_FGS:
                        if checker.check_fg(acyl_fg, product_smiles):
                            if not any(checker.check_fg(acyl_fg, r) for r in reactants_smiles):
                                if checker.check_ring(final_product_heterocycle, product_smiles):
                                    found_acylation = True
                                    acylation_depth = depth
                                    findings_json["atomic_checks"]["functional_groups"].append(acyl_fg)
                                    break

                if not found_halogenation:
                    for halogenation_reaction in HALOGENATION_REACTIONS:
                        if checker.check_reaction(halogenation_reaction, rsmi):
                            if checker.check_ring(final_product_heterocycle, product_smiles):
                                found_halogenation = True
                                halogenation_depth = depth
                                findings_json["atomic_checks"]["named_reactions"].append(halogenation_reaction)
                                break

                if not found_halogenation:
                    for halogen_fg in HALOGEN_FGS:
                        if checker.check_fg(halogen_fg, product_smiles):
                            if not any(checker.check_fg(halogen_fg, r) for r in reactants_smiles):
                                if checker.check_ring(final_product_heterocycle, product_smiles):
                                    found_halogenation = True
                                    halogenation_depth = depth
                                    findings_json["atomic_checks"]["functional_groups"].append(halogen_fg)
                                    break

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = False
    if found_heterocycle and found_acylation and found_halogenation:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "heterocycle_presence",
                    "acylation_event",
                    "halogenation_event"
                ],
                "description": "The route must contain a specific heterocycle, an acylation event, and a halogenation event."
            }
        })
        if heterocycle_depth != -1:
            if heterocycle_depth > acylation_depth > halogenation_depth and halogenation_depth > 0:
                result = True
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "ordered_targets": [
                            "acylation_event",
                            "halogenation_event"
                        ],
                        "description": "The acylation event must occur in an earlier synthetic step than the halogenation event."
                    }
                })
                findings_json["structural_constraints"].append({
                    "type": "sequence",
                    "details": {
                        "ordered_targets": [
                            "ring_formation",
                            "acylation_event"
                        ],
                        "description": "If the target heterocycle is formed de novo during the synthesis, its formation must occur in an earlier step than the acylation event."
                    }
                })
                if halogenation_depth != len(route.get("children", [])) - 1: # Assuming depth 0 is root, and last step is max depth
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "halogenation_event",
                            "position": "not_last_stage",
                            "description": "The halogenation event cannot be the final step of the synthesis."
                        }
                    })

        elif acylation_depth > halogenation_depth and halogenation_depth > 0:
            result = True
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "ordered_targets": [
                        "acylation_event",
                        "halogenation_event"
                    ],
                    "description": "The acylation event must occur in an earlier synthetic step than the halogenation event."
                }
            })
            if halogenation_depth != len(route.get("children", [])) - 1:
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "halogenation_event",
                        "position": "not_last_stage",
                        "description": "The halogenation event cannot be the final step of the synthesis."
                    }
                })

    return result, findings_json
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


BIARYL_REACTION_TYPES = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
    "Suzuki",
    "Stille reaction_aryl",
    "Stille reaction_aryl OTf",
    "Stille",
    "Negishi",
    "Kumada cross-coupling",
    "Hiyama-Denmark Coupling",
    "Aryllithium cross-coupling",
    "Ullmann condensation",
]

AMIDE_REACTION_TYPES = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthetic strategy characterized by an early-stage biaryl bond formation followed by at least two subsequent amide couplings, with one amide coupling occurring as the final step of the synthesis. The biaryl formation is identified by name (using the list of C-C cross-coupling reactions in `BIARYL_REACTION_TYPES`) or by a structural increase in biaryl substructures. The amide couplings are identified by name (using the list in `AMIDE_REACTION_TYPES`) or by a structural increase in amide bonds.
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

    # Track key features
    has_biaryl_forming_reaction = False
    has_biaryl_formation = False
    biaryl_depth = None
    amide_coupling_count = 0
    amide_depths = []
    final_step_is_amide = False

    def dfs_traverse(node, depth=0):
        nonlocal has_biaryl_forming_reaction, has_biaryl_formation, biaryl_depth, amide_coupling_count, amide_depths, final_step_is_amide, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            is_biaryl_reaction = False
            for rxn_type in BIARYL_REACTION_TYPES:
                if checker.check_reaction(rxn_type, rsmi):
                    is_biaryl_reaction = True
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    break

            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol:
                biaryl_patterns = [
                    Chem.MolFromSmarts("c:c-c:c"),  # aromatic-aromatic
                ]

                has_biaryl_in_product = any(
                    product_mol.HasSubstructMatch(pattern) for pattern in biaryl_patterns
                )

                if has_biaryl_in_product:
                    biaryl_in_reactants = False
                    for r in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol and any(
                            r_mol.HasSubstructMatch(pattern) for pattern in biaryl_patterns
                        ):
                            biaryl_in_reactants = True
                            break

                    if not biaryl_in_reactants:
                        has_biaryl_formation = True
                        biaryl_depth = depth
                        findings_json["atomic_checks"]["named_reactions"].append("biaryl_formation")
                        if is_biaryl_reaction:
                            has_biaryl_forming_reaction = True

            is_amide_coupling = False
            for rxn_type in AMIDE_REACTION_TYPES:
                if checker.check_reaction(rxn_type, rsmi):
                    is_amide_coupling = True
                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    break

            if not is_amide_coupling:
                if product_mol:
                    fg_found = []
                    if checker.check_fg("Primary amide", product_smiles):
                        fg_found.append("Primary amide")
                    if checker.check_fg("Secondary amide", product_smiles):
                        fg_found.append("Secondary amide")
                    if checker.check_fg("Tertiary amide", product_smiles):
                        fg_found.append("Tertiary amide")
                    
                    for fg_name in fg_found:
                        if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(fg_name)

                    amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
                    amide_count_in_product = len(product_mol.GetSubstructMatches(amide_pattern))
                    amide_count_in_reactants = 0
                    for r in reactants_smiles:
                        r_mol = Chem.MolFromSmiles(r)
                        if r_mol:
                            amide_count_in_reactants += len(
                                r_mol.GetSubstructMatches(amide_pattern)
                            )
                    if amide_count_in_product > amide_count_in_reactants:
                        is_amide_coupling = True
                        findings_json["atomic_checks"]["named_reactions"].append("amide_formation")

            if is_amide_coupling:
                amide_coupling_count += 1
                amide_depths.append(depth)
                if depth == 1:
                    final_step_is_amide = True

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types where depth increases
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    strategy_detected = (
        (has_biaryl_forming_reaction or has_biaryl_formation)
        and amide_coupling_count >= 2
        and final_step_is_amide
    )

    if strategy_detected and biaryl_depth is not None and amide_depths:
        earliest_amide_depth = min(amide_depths)
        if biaryl_depth <= earliest_amide_depth:
            strategy_detected = False

    # Populate structural constraints based on detected flags
    if has_biaryl_forming_reaction or has_biaryl_formation:
        # This represents the 'biaryl_formation' part of the co-occurrence
        # We don't have a direct structural constraint object for just 'biaryl_formation' in the input JSON,
        # but it's implied by the overall strategy.
        pass # Handled by the overall strategy_detected check

    if amide_coupling_count >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "amide_coupling",
                "operator": ">=",
                "value": 2
            }
        })

    if final_step_is_amide:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "amide_coupling",
                "position": "last_stage"
            }
        })

    if (has_biaryl_forming_reaction or has_biaryl_formation) and amide_coupling_count > 0:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "biaryl_formation",
                    "amide_coupling"
                ]
            }
        })

    if strategy_detected and biaryl_depth is not None and amide_depths:
        earliest_amide_depth = min(amide_depths)
        if biaryl_depth < earliest_amide_depth: # Strictly before
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "biaryl_formation",
                    "after": "amide_coupling"
                }
            })

    return strategy_detected, findings_json
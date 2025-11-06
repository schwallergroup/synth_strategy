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
    "thiazole", "oxazole", "imidazole", "pyrazole", "isoxazole", "isothiazole",
    "triazole", "tetrazole", "pyridine", "pyrimidine", "pyrazine", "pyridazine",
    "furan", "thiophene", "pyrrole", "benzothiazole", "benzoxazole",
    "benzimidazole", "indole", "quinoline", "isoquinoline",
]

RING_FORMATION_REACTIONS = [
    "Formation of NOS Heterocycles", "benzothiazole", "benzoxazole",
    "benzimidazole", "thiazole", "oxazole", "imidazole", "Paal-Knorr pyrrole",
    "Fischer indole", "pyrazole", "indole", "benzofuran", "benzothiophene",
    "oxadiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects syntheses that feature both a late-stage convergent reaction and a late-stage ring-forming reaction. The convergent step is defined as a reaction between two or more complex fragments. The ring-forming step is identified by specific named reactions (e.g., Fischer indole, Paal-Knorr pyrrole), a general increase in ring count, or the formation of a new heterocycle (e.g., pyridine, indole, thiazole). Both events must occur within the final three steps of the synthesis.
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

    # Track if we found the pattern
    found_pattern = False
    # Track if we found a late-stage ring formation
    has_late_ring_formation = False
    # Track if we found a convergent final step
    has_convergent_final = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern, has_late_ring_formation, has_convergent_final, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            # Check if this is a late-stage step (depth 0, 1, or 2)
            if depth <= 2:
                # Record positional constraint for ring formation if applicable
                if not has_late_ring_formation:
                    # Check for ring formation reactions
                    for rxn_type in RING_FORMATION_REACTIONS:
                        if checker.check_reaction(rxn_type, rsmi):
                            has_late_ring_formation = True
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                            if {"type": "positional", "details": {"target": "ring_formation", "position": "within_last_3_stages", "description": "A ring formation event (either a named reaction, an increase in ring count, or de novo heterocycle synthesis) must occur within the last three steps (depth <= 2)."}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ring_formation", "position": "within_last_3_stages", "description": "A ring formation event (either a named reaction, an increase in ring count, or de novo heterocycle synthesis) must occur within the last three steps (depth <= 2)."}})
                            break

                # If no specific reaction type detected, check ring count
                if not has_late_ring_formation:
                    reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r]
                    product_mol = Chem.MolFromSmiles(product_str) if product_str else None

                    if product_mol and all(r for r in reactants_mols):
                        product_rings = product_mol.GetRingInfo().NumRings()
                        reactant_rings_total = sum(r.GetRingInfo().NumRings() for r in reactants_mols)

                        if product_rings > reactant_rings_total:
                            has_late_ring_formation = True
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation") # Generic ring formation
                            if {"type": "positional", "details": {"target": "ring_formation", "position": "within_last_3_stages", "description": "A ring formation event (either a named reaction, an increase in ring count, or de novo heterocycle synthesis) must occur within the last three steps (depth <= 2)."}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ring_formation", "position": "within_last_3_stages", "description": "A ring formation event (either a named reaction, an increase in ring count, or de novo heterocycle synthesis) must occur within the last three steps (depth <= 2)."}})
                        else:
                            for ring in HETEROCYCLES_OF_INTEREST:
                                if checker.check_ring(ring, product_str) and not any(
                                    checker.check_ring(ring, r) for r in reactants_str.split(".")
                                ):
                                    has_late_ring_formation = True
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                                    if {"type": "positional", "details": {"target": "ring_formation", "position": "within_last_3_stages", "description": "A ring formation event (either a named reaction, an increase in ring count, or de novo heterocycle synthesis) must occur within the last three steps (depth <= 2)."}} not in findings_json["structural_constraints"]:
                                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ring_formation", "position": "within_last_3_stages", "description": "A ring formation event (either a named reaction, an increase in ring count, or de novo heterocycle synthesis) must occur within the last three steps (depth <= 2)."}})
                                    break

                # Check if convergent (multiple complex reactants)
                reactants_list = reactants_str.split(".")
                if len(reactants_list) >= 2:
                    complex_reactants = 0
                    for r_smi in reactants_list:
                        r_mol = Chem.MolFromSmiles(r_smi)
                        if r_mol and r_mol.GetNumAtoms() > 6:
                            is_complex = False
                            if r_mol.GetRingInfo().NumRings() > 0:
                                is_complex = True
                            else:
                                for fg in [
                                    "Carboxylic acid", "Ester", "Amide", "Amine",
                                    "Alcohol", "Aldehyde", "Ketone",
                                ]:
                                    if checker.check_fg(fg, r_smi):
                                        findings_json["atomic_checks"]["functional_groups"].append(fg)
                                        is_complex = True
                                        break
                            if is_complex:
                                complex_reactants += 1

                    if complex_reactants >= 2:
                        has_convergent_final = True
                        if {"type": "count", "details": {"target": "complex_reactants_per_step", "operator": ">=", "value": 2, "description": "A convergent reaction is defined as a step with at least two complex reactants. A reactant is complex if it has >6 atoms and contains a ring or a key functional group."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "complex_reactants_per_step", "operator": ">=", "value": 2, "description": "A convergent reaction is defined as a step with at least two complex reactants. A reactant is complex if it has >6 atoms and contains a ring or a key functional group."}})
                        if {"type": "positional", "details": {"target": "convergent_reaction", "position": "within_last_3_stages", "description": "A convergent reaction, defined as a step with at least two complex reactants, must occur within the last three steps (depth <= 2)."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "convergent_reaction", "position": "within_last_3_stages", "description": "A convergent reaction, defined as a step with at least two complex reactants, must occur within the last three steps (depth <= 2)."}})

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if both required elements of the pattern are present
    found_pattern = has_late_ring_formation and has_convergent_final

    if found_pattern:
        if {"type": "co-occurrence", "details": {"targets": ["late_stage_ring_formation", "late_stage_convergent_reaction"], "description": "The synthesis must contain both a ring formation event and a convergent reaction event, with both occurring in the final three steps."}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["late_stage_ring_formation", "late_stage_convergent_reaction"], "description": "The synthesis must contain both a ring formation event and a convergent reaction event, with both occurring in the final three steps."}})

    # Deduplicate lists in findings_json
    for key in findings_json["atomic_checks"]:
        findings_json["atomic_checks"][key] = sorted(list(set(findings_json["atomic_checks"][key])))
    findings_json["structural_constraints"] = sorted(list(set(tuple(sorted(d.items())) for d in findings_json["structural_constraints"])), key=lambda x: str(x))
    findings_json["structural_constraints"] = [dict(t) for t in findings_json["structural_constraints"]]

    return found_pattern, findings_json

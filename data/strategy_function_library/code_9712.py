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


HETEROAROMATIC_RINGS_OF_INTEREST = [
    "pyrazole",
    "triazole",
    "tetrazole",
    "imidazole",
    "oxazole",
    "thiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for synthetic routes that feature the co-occurrence of several key events: 1) the use of a hydrazine-containing reactant, 2) the use of a carbonyl-containing reactant (Aldehyde, Ketone, Ester, or Carboxylic acid), 3) the formation of a specific heteroaromatic ring (pyrazole, triazole, tetrazole, imidazole, oxazole, or thiazole), and 4) a late-stage decarboxylation occurring in the final two steps of the synthesis.
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

    # Track if we found the key features
    found_hydrazine = False
    found_carbonyl = False
    found_ring_formation = False
    found_aromatization = False
    found_decarboxylation = False
    decarboxylation_depth = float("inf")
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal found_hydrazine, found_carbonyl, found_ring_formation, found_aromatization
        nonlocal found_decarboxylation, decarboxylation_depth, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for hydrazine in reactants
            for r_smiles in reactants_smiles:
                if checker.check_fg("Hydrazine", r_smiles):
                    found_hydrazine = True
                    if "Hydrazine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Hydrazine")

            # Check for carbonyl compounds in reactants
            for r_smiles in reactants_smiles:
                if checker.check_fg("Aldehyde", r_smiles):
                    found_carbonyl = True
                    if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                if checker.check_fg("Ketone", r_smiles):
                    found_carbonyl = True
                    if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                if checker.check_fg("Ester", r_smiles):
                    found_carbonyl = True
                    if "Ester" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Ester")
                if checker.check_fg("Carboxylic acid", r_smiles):
                    found_carbonyl = True
                    if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")

            # Specifically check for the formation of certain heteroaromatic rings
            for ring in HETEROAROMATIC_RINGS_OF_INTEREST:
                if checker.check_ring(ring, product_smiles) and not any(checker.check_ring(ring, r) for r in reactants_smiles):
                    found_ring_formation = True
                    found_aromatization = True # All rings in the list are aromatic
                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                    if "aromatization" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("aromatization")

            # Check for decarboxylation
            is_decarboxylation_reaction = checker.check_reaction("Decarboxylation", rsmi)
            is_decarboxylation_fg_change = (
                any(checker.check_fg("Carboxylic acid", r) for r in reactants_smiles)
                and not checker.check_fg("Carboxylic acid", product_smiles)
            )

            if is_decarboxylation_reaction or is_decarboxylation_fg_change:
                found_decarboxylation = True
                decarboxylation_depth = min(decarboxylation_depth, depth)
                if "Decarboxylation" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Decarboxylation")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for chemical children
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for reaction children
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    is_late_stage_decarboxylation = (decarboxylation_depth <= 2)

    strategy_present = (
        found_hydrazine
        and found_carbonyl
        and found_ring_formation
        and found_aromatization
        and found_decarboxylation
        and is_late_stage_decarboxylation
    )

    # Populate structural constraints based on the final flags
    if found_hydrazine and found_carbonyl and found_ring_formation and found_aromatization and found_decarboxylation:
        # This corresponds to the first structural constraint (co-occurrence)
        if {"type": "co-occurrence", "details": {"targets": ["Hydrazine", "carbonyl_group_reactant", "heteroaromatic_ring_formation", "aromatization", "Decarboxylation"]}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Hydrazine", "carbonyl_group_reactant", "heteroaromatic_ring_formation", "aromatization", "Decarboxylation"]}})

    if found_decarboxylation and is_late_stage_decarboxylation:
        # This corresponds to the second structural constraint (positional)
        if {"type": "positional", "details": {"target": "Decarboxylation", "position": "last_two_steps"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Decarboxylation", "position": "last_two_steps"}})

    return strategy_present, findings_json

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


CYCLIC_ACETAL_KETAL_RINGS = ["dioxolane", "dioxane", "dioxolene", "dioxepane", "trioxane"]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route involves protection of diols as cyclic acetals/ketals.
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

    # Track if we found a valid protection reaction
    found_protection_reaction = False

    def dfs_traverse(node):
        nonlocal found_protection_reaction, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for protection reactions (acetalization)
            is_diol_acetalization = checker.check_reaction("Diol acetalization", rsmi)
            is_aldehyde_ketone_acetalization = checker.check_reaction("Aldehyde or ketone acetalization", rsmi)

            if is_diol_acetalization or is_aldehyde_ketone_acetalization:
                if is_diol_acetalization:
                    if "Diol acetalization" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Diol acetalization")
                if is_aldehyde_ketone_acetalization:
                    if "Aldehyde or ketone acetalization" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Aldehyde or ketone acetalization")

                has_cyclic_acetal = False
                for ring in CYCLIC_ACETAL_KETAL_RINGS:
                    if checker.check_ring(ring, product):
                        has_cyclic_acetal = True
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)

                if has_cyclic_acetal:
                    found_protection_reaction = True
                    # Add structural constraints for acetalization reactions
                    if is_diol_acetalization:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Diol acetalization", "dioxolane"]}})
                    if is_aldehyde_ketone_acetalization:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Aldehyde or ketone acetalization", "dioxolane"]}})
                    # No return here, continue to check other conditions if needed for more findings

            # Direct detection of cyclic acetal/ketal formation
            has_cyclic_acetal_product = False
            for ring in CYCLIC_ACETAL_KETAL_RINGS:
                if checker.check_ring(ring, product):
                    has_cyclic_acetal_product = True
                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(ring)

            has_cyclic_acetal_reactants = False
            for reactant in reactants:
                for ring in CYCLIC_ACETAL_KETAL_RINGS:
                    if checker.check_ring(ring, reactant):
                        has_cyclic_acetal_reactants = True
                        # Do not add to findings_json if it's in reactants, as it's not 'formed'

            # If a cyclic acetal/ketal is formed (appears in product but not in reactants)
            if has_cyclic_acetal_product and not has_cyclic_acetal_reactants:
                has_diol = False
                has_carbonyl = False
                alcohol_count_total = 0

                for reactant in reactants:
                    # Count alcohols in this specific reactant
                    alcohol_count = 0
                    if checker.check_fg("Primary alcohol", reactant):
                        alcohol_count += 1
                        if "Primary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary alcohol")
                    if checker.check_fg("Secondary alcohol", reactant):
                        alcohol_count += 1
                        if "Secondary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary alcohol")
                    if checker.check_fg("Tertiary alcohol", reactant):
                        alcohol_count += 1
                        if "Tertiary alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary alcohol")

                    # If this single reactant has multiple alcohols, it's a potential diol
                    if alcohol_count >= 2:
                        has_diol = True
                        alcohol_count_total += alcohol_count # Accumulate for structural constraint
                        if {"type": "count", "details": {"target": "alcohol_groups_on_single_reactant", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "alcohol_groups_on_single_reactant", "operator": ">=", "value": 2}})

                    # Check for carbonyl compounds
                    if checker.check_fg("Aldehyde", reactant):
                        has_carbonyl = True
                        if "Aldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                    if checker.check_fg("Ketone", reactant):
                        has_carbonyl = True
                        if "Ketone" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Ketone")
                    if checker.check_fg("Formaldehyde", reactant):
                        has_carbonyl = True
                        if "Formaldehyde" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Formaldehyde")

                # If we have either a diol or a carbonyl, and a cyclic acetal is formed, it's likely a protection
                if has_diol or has_carbonyl:
                    found_protection_reaction = True
                    # Add structural constraints for ring formation with diol/carbonyl
                    if has_diol and {"type": "co-occurrence", "details": {"targets": ["ring_formation", "Primary alcohol"]}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["ring_formation", "Primary alcohol"]}})
                    if has_carbonyl and {"type": "co-occurrence", "details": {"targets": ["ring_formation", "Aldehyde"]}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["ring_formation", "Aldehyde"]}})
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

            # Check for deprotection reactions (hydrolysis)
            is_acetal_hydrolysis_diol = checker.check_reaction("Acetal hydrolysis to diol", rsmi)
            is_acetal_hydrolysis_aldehyde = checker.check_reaction("Acetal hydrolysis to aldehyde", rsmi)
            is_ketal_hydrolysis_ketone = checker.check_reaction("Ketal hydrolysis to ketone", rsmi)

            if is_acetal_hydrolysis_diol or is_acetal_hydrolysis_aldehyde or is_ketal_hydrolysis_ketone:
                if is_acetal_hydrolysis_diol:
                    if "Acetal hydrolysis to diol" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Acetal hydrolysis to diol")
                if is_acetal_hydrolysis_aldehyde:
                    if "Acetal hydrolysis to aldehyde" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Acetal hydrolysis to aldehyde")
                if is_ketal_hydrolysis_ketone:
                    if "Ketal hydrolysis to ketone" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Ketal hydrolysis to ketone")

                has_cyclic_acetal_reactant = False
                for reactant in reactants:
                    for ring in CYCLIC_ACETAL_KETAL_RINGS:
                        if checker.check_ring(ring, reactant):
                            has_cyclic_acetal_reactant = True
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)

                if has_cyclic_acetal_reactant:
                    found_protection_reaction = True
                    # Add structural constraints for hydrolysis reactions
                    if is_acetal_hydrolysis_diol:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Acetal hydrolysis to diol", "dioxolane"]}})
                    if is_acetal_hydrolysis_aldehyde:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Acetal hydrolysis to aldehyde", "dioxolane"]}})
                    if is_ketal_hydrolysis_ketone:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Ketal hydrolysis to ketone", "dioxolane"]}})

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return found_protection_reaction, findings_json

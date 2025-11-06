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
    "pyridine",
    "furan",
    "thiophene",
    "pyrrole",
    "imidazole",
    "oxazole",
    "thiazole",
    "pyrazole",
    "isoxazole",
    "isothiazole",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "triazole",
    "tetrazole",
]

ARYLATION_REACTIONS_OF_INTEREST = [
    "Suzuki",
    "Negishi",
    "Buchwald-Hartwig",
    "N-arylation",
    "Heck",
    "Stille",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Ullmann-Goldberg Substitution amine",
    "Ullmann-Goldberg Substitution thiol",
    "Ullmann-Goldberg Substitution aryl alcohol",
    "Ullmann condensation",
    "Chan-Lam amine",
    "Chan-Lam alcohol",
    "Chan-Lam etherification",
]

ARYL_GROUPS_OF_INTEREST = ["Phenol", "Aniline", "Aromatic halide"]

def has_aryl_group(smiles, findings_json_ref):
    """Check if a molecule contains an aryl group"""
    for aryl in ARYL_GROUPS_OF_INTEREST:
        if checker.check_fg(aryl, smiles):
            print(f"Found aryl group: {aryl} in {smiles}")
            if aryl not in findings_json_ref["atomic_checks"]["functional_groups"]:
                findings_json_ref["atomic_checks"]["functional_groups"].append(aryl)
            return True
    # Check for benzene ring as a fallback
    if checker.check_ring("benzene", smiles):
        print(f"Found benzene ring in {smiles}")
        if "benzene" not in findings_json_ref["atomic_checks"]["ring_systems"]:
            findings_json_ref["atomic_checks"]["ring_systems"].append("benzene")
        return True
    return False

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage arylation strategy. The function specifically checks for reactions
    occurring in the final two steps of a synthesis that match a predefined list of
    arylation reaction types (e.g., Suzuki, Buchwald-Hartwig) defined in
    `ARYLATION_REACTIONS_OF_INTEREST`. It confirms that one of the reactants is a
    heterocycle from the `HETEROCYCLES_OF_INTEREST` list and that this heterocycle
    is preserved in the product.
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

    found_late_arylation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_arylation, findings_json

        if node["type"] == "reaction" and depth <= 2:  # Check reactions up to depth 2 (late stage)
            try:
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    return

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is an arylation reaction
                is_arylation = False
                for reaction_type in ARYLATION_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_arylation = True
                        print(f"Found arylation reaction type: {reaction_type}")
                        if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

                if not is_arylation:
                    print("Not an arylation reaction")
                    return
                
                # Structural constraint: positional (depth <= 2)
                if {"type": "positional", "details": {"target": "arylation_reaction", "position": "depth <= 2"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "arylation_reaction", "position": "depth <= 2"}})

                # Check for aryl groups in reactants
                has_aryl_reactant = False
                for reactant in reactants_smiles:
                    if has_aryl_group(reactant, findings_json):
                        has_aryl_reactant = True
                        print(f"Found aryl group in reactant: {reactant}")
                        break

                if not has_aryl_reactant:
                    print("No aryl group found in reactants")
                    return

                # Check for heterocycles in reactants
                heterocycle_reactant = None
                heterocycle_found = None
                for reactant in reactants_smiles:
                    for heterocycle in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(heterocycle, reactant):
                            heterocycle_reactant = reactant
                            heterocycle_found = heterocycle
                            print(f"Found heterocycle {heterocycle} in reactant: {reactant}")
                            if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                            break
                    if heterocycle_reactant:
                        break

                # Check if product contains both the heterocycle and aryl group
                if heterocycle_found:
                    if checker.check_ring(heterocycle_found, product_smiles):
                        print(f"Product contains heterocycle {heterocycle_found}")
                        if heterocycle_found not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(heterocycle_found)
                        print(f"Confirmed late-stage heterocycle arylation at depth {depth}")
                        found_late_arylation = True
                        
                        # Structural constraint: co-occurrence
                        co_occurrence_constraint = {"type": "co-occurrence", "details": {"targets": ["arylation_reaction", "reactant_with_heterocycle", "reactant_with_aryl_group"], "scope": "single_reaction_step"}}
                        if co_occurrence_constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(co_occurrence_constraint)

                        # Structural constraint: sequence
                        sequence_constraint = {"type": "sequence", "details": {"before": "heterocycle_in_reactant", "after": "same_heterocycle_in_product", "scope": "single_reaction_step"}}
                        if sequence_constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(sequence_constraint)

                    else:
                        print(f"Product does not contain the heterocycle {heterocycle_found}")
                else:
                    print("No heterocycle found in reactants")
            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        # Process children with increased depth
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {found_late_arylation}")
    return found_late_arylation, findings_json

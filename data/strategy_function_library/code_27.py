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
    This function detects a strategy involving late-stage nitrogen functionalization.
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

    # Track if we found the key reactions
    n_functionalization_steps = []
    total_steps = 0

    print("Starting analysis for late-stage N-functionalization strategy")

    def dfs_traverse(node, depth=0):
        nonlocal n_functionalization_steps, total_steps, findings_json

        if node["type"] == "reaction":
            total_steps += 1
            print(f"Processing reaction at depth {depth}")

            # Extract reactants and product
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Reaction SMILES: {rsmi}")

                # Check for N-functionalization reactions
                is_n_functionalization = False

                # First verify that the product contains nitrogen
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol:
                    has_nitrogen = any(atom.GetAtomicNum() == 7 for atom in product_mol.GetAtoms())

                    if has_nitrogen:
                        print(f"Product contains nitrogen at depth {depth}")

                        # Check for N-alkylation reactions
                        n_alkylation_reactions = [
                            "N-alkylation of primary amines with alkyl halides",
                            "N-alkylation of secondary amines with alkyl halides",
                            "Methylation with MeI_primary",
                            "Methylation with MeI_secondary",
                            "Methylation with MeI_tertiary",
                            "DMS Amine methylation",
                            "Eschweiler-Clarke Primary Amine Methylation",
                            "Eschweiler-Clarke Secondary Amine Methylation",
                            "Reductive methylation of primary amine with formaldehyde",
                        ]
                        for r_name in n_alkylation_reactions:
                            if checker.check_reaction(r_name, rsmi):
                                is_n_functionalization = True
                                findings_json["atomic_checks"]["named_reactions"].append(r_name)
                                print(f"Found N-alkylation: {r_name} at depth {depth}")
                                # Additional check for methylation specifically
                                if "I[CH3" in rsmi or "[CH3:10]" in rsmi:
                                    print(f"Found methylation pattern in SMILES at depth {depth}")
                                break

                        # Check for N-acylation reactions
                        if not is_n_functionalization:
                            n_acylation_reactions = [
                                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                                "Acyl chloride with secondary amine to amide",
                                "Carboxylic acid with primary amine to amide",
                                "Ester with primary amine to amide",
                                "Ester with secondary amine to amide",
                                "{Schotten-Baumann_amide}",
                                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                            ]
                            for r_name in n_acylation_reactions:
                                if checker.check_reaction(r_name, rsmi):
                                    is_n_functionalization = True
                                    findings_json["atomic_checks"]["named_reactions"].append(r_name)
                                    print(f"Found N-acylation: {r_name} at depth {depth}")
                                    break

                        # Check for N-arylation reactions
                        if not is_n_functionalization:
                            n_arylation_reactions = [
                                "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                                "{Buchwald-Hartwig}",
                                "{N-arylation_heterocycles}",
                                "Goldberg coupling",
                                "Ullmann-Goldberg Substitution amine",
                            ]
                            for r_name in n_arylation_reactions:
                                if checker.check_reaction(r_name, rsmi):
                                    is_n_functionalization = True
                                    findings_json["atomic_checks"]["named_reactions"].append(r_name)
                                    print(f"Found N-arylation: {r_name} at depth {depth}")
                                    break

                        # Check for sulfonamide formation
                        if not is_n_functionalization:
                            sulfonamide_reactions = [
                                "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                                "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
                                "{sulfon_amide}",
                            ]
                            for r_name in sulfonamide_reactions:
                                if checker.check_reaction(r_name, rsmi):
                                    is_n_functionalization = True
                                    findings_json["atomic_checks"]["named_reactions"].append(r_name)
                                    print(f"Found sulfonamide formation: {r_name} at depth {depth}")
                                    break

                        # Check for urea/thiourea formation
                        if not is_n_functionalization:
                            urea_thiourea_reactions = [
                                "Urea synthesis via isocyanate and primary amine",
                                "Urea synthesis via isocyanate and secondary amine",
                                "{urea}",
                                "{thiourea}",
                            ]
                            for r_name in urea_thiourea_reactions:
                                if checker.check_reaction(r_name, rsmi):
                                    is_n_functionalization = True
                                    findings_json["atomic_checks"]["named_reactions"].append(r_name)
                                    print(f"Found urea/thiourea formation: {r_name} at depth {depth}")
                                    break

                        # Check for reductive amination
                        if not is_n_functionalization:
                            reductive_amination_reactions = [
                                "Reductive amination with aldehyde",
                                "Reductive amination with ketone",
                                "Reductive amination with alcohol",
                                "{reductive amination}",
                            ]
                            for r_name in reductive_amination_reactions:
                                if checker.check_reaction(r_name, rsmi):
                                    is_n_functionalization = True
                                    findings_json["atomic_checks"]["named_reactions"].append(r_name)
                                    print(f"Found reductive amination: {r_name} at depth {depth}")
                                    break

                        # If we found a N-functionalization reaction, record its depth
                        if is_n_functionalization:
                            n_functionalization_steps.append(depth)
                            print(f"Confirmed N-functionalization at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            # Depth remains the same when traversing from reaction to chemical node
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Total steps: {total_steps}")
    print(f"N-functionalization steps at depths: {n_functionalization_steps}")

    # Check if N-functionalization occurs in the late stage of the synthesis
    # In retrosynthesis, lower depths correspond to later stages in forward synthesis
    result = False
    if total_steps > 0 and n_functionalization_steps:
        # Consider the first half of steps as late stage (low depth in retrosynthesis)
        late_stage_threshold = total_steps // 2

        # In retrosynthesis, late stage reactions have low depth values
        late_stage_n_functionalizations = sum(
            1 for d in n_functionalization_steps if d <= late_stage_threshold
        )

        print(f"Late stage threshold: {late_stage_threshold}")
        print(f"Late stage N-functionalizations: {late_stage_n_functionalizations}")

        if late_stage_n_functionalizations > 0:
            result = True
            # Add the structural constraint if the condition is met
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "any_N_functionalization_reaction",
                    "position": "late_stage",
                    "definition": "At least one N-functionalization reaction must occur within the first half of the retrosynthetic steps (depth <= total_steps / 2)."
                }
            })

    return result, findings_json

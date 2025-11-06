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


SNAR_ACTIVATING_HETEROCYCLES = [
    "pyridine",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "triazole",
    "tetrazole",
    "oxazole",
    "thiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route involves a sequential strategy of at least two
    nucleophilic aromatic substitution (SNAr) or related C-X cross-coupling
    reactions. This includes named reactions like Buchwald-Hartwig and Ullmann,
    and a general pattern for SNAr on rings activated by electron-withdrawing
    groups or specific heterocycles (e.g., pyridine, pyrimidine, pyrazine,
    pyridazine, triazole, tetrazole, oxazole, thiazole).
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

    snar_reactions = []
    result = False

    def calculate_depth(current_node, root_node, current_depth=0):
        """Helper function to calculate depth if not provided in metadata"""
        if current_node == root_node:
            return 0

        # Traverse the tree to find the current node's depth
        def find_node_depth(node, target, depth=0):
            if node == target:
                return depth

            for child in node.get("children", []):
                result = find_node_depth(child, target, depth + 1)
                if result is not None:
                    return result
            return None

        return find_node_depth(root_node, current_node)

    def dfs_traverse(node, current_depth=0):
        nonlocal result
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Get depth from metadata or calculate it
            if "depth" in node["metadata"]:
                depth = int(node["metadata"]["depth"])
            else:
                depth = current_depth

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check for nucleophilic aromatic substitution reactions using predefined reactions
            is_snar = False
            snar_reaction_names = [
                "heteroaromatic_nuc_sub",
                "nucl_sub_aromatic_ortho_nitro",
                "nucl_sub_aromatic_para_nitro",
                "N-arylation",
                "Buchwald-Hartwig",
                "Ullmann-Goldberg Substitution amine",
                "Ullmann-Goldberg Substitution thiol",
                "Ullmann-Goldberg Substitution aryl alcohol",
                "Ullmann condensation",
            ]
            for r_name in snar_reaction_names:
                if checker.check_reaction(r_name, rsmi):
                    is_snar = True
                    findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    break

            if is_snar:
                print(f"Found SNAr reaction through reaction check at depth {depth}")
                snar_reactions.append(depth)
            else:
                # If not found in predefined reactions, check for the pattern manually
                print("Checking for SNAr pattern manually...")

                aromatic_with_leaving_group = False
                nucleophile_present = False
                has_ewg = False

                aromatic_rings = [
                    "benzene", "pyridine", "pyrimidine", "pyrazine", "pyridazine",
                    "furan", "thiophene", "pyrrole", "imidazole", "oxazole",
                    "thiazole", "triazole", "tetrazole", "indole", "quinoline",
                    "isoquinoline",
                ]

                for reactant in reactants:
                    try:
                        for ring in aromatic_rings:
                            if checker.check_ring(ring, reactant):
                                findings_json["atomic_checks"]["ring_systems"].append(ring)
                                has_aromatic_ring = True

                        if has_aromatic_ring:
                            print(f"Found aromatic ring in reactant: {reactant}")

                            leaving_groups = [
                                "Aromatic halide", "Triflate", "Tosylate",
                                "Mesylate", "Thiocyanate"
                            ]
                            for fg_name in leaving_groups:
                                if checker.check_fg(fg_name, reactant):
                                    print(f"Found leaving group in reactant: {reactant}")
                                    findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                                    aromatic_with_leaving_group = True

                            ewg_groups = [
                                "Nitro group", "Nitrile", "Ester", "Ketone",
                                "Trifluoro group", "Sulfone", "Sulfonamide"
                            ]
                            for fg_name in ewg_groups:
                                if checker.check_fg(fg_name, reactant):
                                    print(f"Found electron-withdrawing group in reactant: {reactant}")
                                    findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                                    has_ewg = True

                        nucleophile_groups = [
                            "Primary amine", "Secondary amine", "Phenol",
                            "Primary alcohol", "Secondary alcohol", "Aromatic thiol",
                            "Aliphatic thiol", "Aniline"
                        ]
                        for fg_name in nucleophile_groups:
                            if checker.check_fg(fg_name, reactant):
                                print(f"Found nucleophile in reactant: {reactant}")
                                findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                                nucleophile_present = True
                    except Exception as e:
                        print(f"Error analyzing reactant: {e}")
                        continue

                product_has_substitution = False
                try:
                    for ring in aromatic_rings:
                        if checker.check_ring(ring, product):
                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                            has_product_aromatic_ring = True

                    if has_product_aromatic_ring:
                        print(f"Found aromatic ring in product: {product}")
                        product_substitution_groups = [
                            "Aniline", "Phenol", "Ether", "Monosulfide",
                            "Secondary amine", "Tertiary amine"
                        ]
                        for fg_name in product_substitution_groups:
                            if checker.check_fg(fg_name, product):
                                print(f"Found nucleophilic substitution in product: {product}")
                                findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                                product_has_substitution = True
                except Exception as e:
                    print(f"Error analyzing product: {e}")

                is_heterocycle = False
                for reactant in reactants:
                    for ring in SNAR_ACTIVATING_HETEROCYCLES:
                        if checker.check_ring(ring, reactant):
                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                            is_heterocycle = True
                            break
                    if is_heterocycle: # Break outer loop if heterocycle found in any reactant
                        break

                manual_snar = (
                    aromatic_with_leaving_group
                    and nucleophile_present
                    and product_has_substitution
                    and (has_ewg or is_heterocycle)
                )

                if manual_snar:
                    print(f"Found SNAr reaction through manual check at depth {depth}")
                    snar_reactions.append(depth)

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, current_depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, current_depth + 1)

    print("Starting traversal of synthesis route")
    dfs_traverse(route)

    if len(snar_reactions) >= 2:
        snar_reactions.sort()
        print(f"Found SNAr reactions at depths: {snar_reactions}")

        if len(set(snar_reactions)) >= 2:
            print("Confirmed sequential SNAr strategy")
            result = True
            # Add the structural constraint if the condition is met
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "reaction_stage_with_SNAr_or_coupling",
                    "operator": ">=",
                    "value": 2
                }
            })
        else:
            print("Found multiple SNAr reactions but they are at the same depth")
    else:
        print(
            f"Found only {len(snar_reactions)} SNAr reactions, need at least 2 for a sequential strategy"
        )

    # Ensure unique entries in atomic checks lists
    for key in findings_json["atomic_checks"]:
        findings_json["atomic_checks"][key] = list(set(findings_json["atomic_checks"][key]))

    return result, findings_json

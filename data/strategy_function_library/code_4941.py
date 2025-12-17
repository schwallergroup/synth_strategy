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


from rdkit import Chem
from rdkit.Chem import Descriptors

# Assume 'checker' is a pre-defined object with methods like check_reaction, check_fg, etc.

ETHER_FORMATION_REACTIONS = [
    "Williamson Ether Synthesis",
    "Williamson Ether Synthesis (intra to epoxy)",
    "Mitsunobu aryl ether",
    "Chan-Lam etherification",
    "{Williamson ether}",
    "Alcohol to ether",
    "Ullmann-Goldberg Substitution aryl alcohol",
]

COMPLEXITY_HETEROCYCLES = [
    "furan", "pyrrole", "pyridine", "pyrazole", "imidazole", "oxazole",
    "thiazole", "pyrimidine", "piperidine", "morpholine", "thiophene",
    "indole", "benzoxazole", "benzothiazole", "benzimidazole", "quinoline",
    "isoquinoline", "pyrazine", "pyridazine", "triazole", "tetrazole",
]

COMPLEXITY_FGS = [
    "Carboxylic acid", "Ester", "Amide", "Nitrile", "Nitro group",
    "Sulfonamide", "Urea", "Carbamic ester", "Halogen", "Amine",
    "Alcohol", "Ketone", "Aldehyde", "Sulfone", "Sulfoxide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis where two complex fragments are joined via an ether linkage
    in the final step of the synthesis.
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

    # Original strategy JSON for structural constraints reference
    strategy_json_constraints = [
        {
            "type": "positional",
            "details": {
                "target": "ether_formation",
                "position": "last_stage"
            }
        },
        {
            "type": "count",
            "details": {
                "target": "complex_reactant",
                "operator": ">=",
                "value": 2,
                "scope": "last_stage_reactants"
            }
        },
        {
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "last_stage_ether_formation",
                    "two_or_more_complex_reactants"
                ],
                "scope": "last_stage_reaction"
            }
        }
    ]

    final_step_joins_complex_fragments = False
    has_complex_fragment1 = False
    has_complex_fragment2 = False
    has_late_stage_ether_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_joins_complex_fragments, has_complex_fragment1, has_complex_fragment2, has_late_stage_ether_formation, findings_json

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and depth == 1:  # First reaction (late stage)
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing late-stage reaction: {rsmi}")
                print(f"Number of reactants: {len(reactants)}")

                # Check if we have at least 2 reactants
                if len(reactants) >= 2:
                    # Check for ether formation reaction
                    is_ether_formation = False
                    for reaction_name in ETHER_FORMATION_REACTIONS:
                        if checker.check_reaction(reaction_name, rsmi):
                            is_ether_formation = True
                            has_late_stage_ether_formation = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                            break

                    if is_ether_formation:
                        print(f"Ether formation reaction detected: {is_ether_formation}")
                        # Add positional constraint for ether formation
                        findings_json["structural_constraints"].append(strategy_json_constraints[0])
                    else:
                        # Check if product has ether that wasn't in reactants
                        try:
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol and checker.check_fg("Ether", product):
                                findings_json["atomic_checks"]["functional_groups"].append("Ether")
                                # Count ethers in reactants and product
                                product_ether_count = (
                                    len(checker.get_fg_atom_indices("Ether", product))
                                    if checker.check_fg("Ether", product)
                                    else 0
                                )
                                reactant_ether_count = sum(
                                    (
                                        len(checker.get_fg_atom_indices("Ether", r))
                                        if checker.check_fg("Ether", r)
                                        else 0
                                    )
                                    for r in reactants
                                )

                                print(f"Product ether count: {product_ether_count}")
                                print(f"Reactant ether count: {reactant_ether_count}")

                                if product_ether_count > reactant_ether_count:
                                    has_late_stage_ether_formation = True
                                    findings_json["atomic_checks"]["named_reactions"].append("ether_formation") # Generic ether formation
                                    findings_json["structural_constraints"].append(strategy_json_constraints[0])
                                    print("New ether detected based on count")
                                    print("New ether bond detected")
                        except Exception as e:
                            print(f"Error checking ether formation: {e}")

                    # Check complexity of reactants
                    complex_fragments = 0
                    complex_reactants = []

                    for r in reactants:
                        try:
                            r_mol = Chem.MolFromSmiles(r)
                            if not r_mol:
                                continue

                            # Calculate complexity based on multiple factors
                            is_complex = False

                            # 1. Check for heterocycles
                            heterocycle_present = False
                            for ring in COMPLEXITY_HETEROCYCLES:
                                if checker.check_ring(ring, r):
                                    heterocycle_present = True
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                            
                            # 2. Check for ring count
                            ring_info = r_mol.GetRingInfo()
                            ring_count = ring_info.NumRings()

                            # 3. Check for functional groups that add complexity
                            fg_count = 0
                            for fg in COMPLEXITY_FGS:
                                if checker.check_fg(fg, r):
                                    fg_count += 1
                                    findings_json["atomic_checks"]["functional_groups"].append(fg)

                            # 4. Check molecular weight as a proxy for complexity
                            mol_weight = Descriptors.MolWt(r_mol)

                            # 5. Check atom count
                            atom_count = r_mol.GetNumAtoms()

                            # Determine if complex based on combined factors - more lenient criteria
                            is_complex = (
                                heterocycle_present
                                or ring_count >= 1
                                or fg_count >= 2
                                or mol_weight > 120
                                or atom_count > 10
                            )

                            if is_complex:
                                complex_fragments += 1
                                complex_reactants.append(r)
                                print(f"Complex fragment detected: {r}")
                                print(f"  - Has heterocycle: {heterocycle_present}")
                                print(f"  - Ring count: {ring_count}")
                                print(f"  - Functional group count: {fg_count}")
                                print(f"  - Molecular weight: {mol_weight}")
                                print(f"  - Atom count: {atom_count}")
                        except Exception as e:
                            print(f"Error checking fragment complexity: {e}")

                    has_complex_fragment1 = complex_fragments >= 1
                    has_complex_fragment2 = complex_fragments >= 2

                    if has_complex_fragment2:
                        # Add count constraint for complex reactants
                        findings_json["structural_constraints"].append(strategy_json_constraints[1])

                    # If we have at least two complex fragments and ether formation,
                    # mark as convergent synthesis with late-stage ether linkage
                    if has_complex_fragment2 and has_late_stage_ether_formation:
                        print(
                            "Found at least two complex fragments and ether formation in final step"
                        )
                        final_step_joins_complex_fragments = True
                        # Add co-occurrence constraint
                        findings_json["structural_constraints"].append(strategy_json_constraints[2])

        # Process children
        for child in node.get("children", []):
            new_depth = depth + 1 if node['type'] != 'reaction' else depth
            dfs_traverse(child, new_depth)

    # Start traversal from the root node (final product)
    if route["type"] == "mol":
        print(f"Starting traversal from final product: {route['smiles']}")
        dfs_traverse(route)
    else:
        print("Warning: Expected root node to be a molecule")

    # Final result
    print(f"Convergent synthesis detection: {final_step_joins_complex_fragments}")
    print(f"  - Has complex fragment 1: {has_complex_fragment1}")
    print(f"  - Has complex fragment 2: {has_complex_fragment2}")
    print(f"  - Has late-stage ether formation: {has_late_stage_ether_formation}")

    # Remove duplicates from atomic checks lists
    for key in findings_json["atomic_checks"]:
        findings_json["atomic_checks"][key] = list(set(findings_json["atomic_checks"][key]))
    
    # Remove duplicate structural constraints (based on content)
    unique_constraints = []
    seen_constraints = set()
    for constraint in findings_json["structural_constraints"]:
        constraint_str = str(constraint) # Convert dict to string for set comparison
        if constraint_str not in seen_constraints:
            unique_constraints.append(constraint)
            seen_constraints.add(constraint_str)
    findings_json["structural_constraints"] = unique_constraints

    return final_step_joins_complex_fragments, findings_json

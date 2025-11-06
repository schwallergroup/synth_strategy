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


AROMATIC_RINGS_OF_INTEREST = [
    "benzene",
    "naphthalene",
    "anthracene",
    "pyridine",
    "furan",
    "thiophene",
    "pyrrole",
    "indole",
    "quinoline",
    "isoquinoline",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis strategy that connects specific aromatic fragments via ether linkages and uses a late-stage esterification. The strategy is identified by finding at least one convergent step, a late-stage (depth <= 3) esterification reaction, and one or more ether-forming reactions that involve an aromatic reactant from a predefined list (e.g., benzene, pyridine, furan). It also confirms the presence of a polyether structure (>=2 ether groups) in reaction products.
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

    # Initialize tracking variables
    ether_formations = 0
    late_esterification = False
    has_polyether_linker = False
    has_aromatic_fragments = 0
    branching_points = 0
    final_product_has_ester = False

    # Check if final product has ester group
    if route["type"] == "mol" and checker.check_fg("Ester", route["smiles"]):
        final_product_has_ester = True
        findings_json["atomic_checks"]["functional_groups"].append("Ester")
        print(f"Final product has ester functional group: {route['smiles']}")

    def dfs_traverse(node, depth=0):
        nonlocal ether_formations, late_esterification, has_polyether_linker, has_aromatic_fragments, branching_points, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for late-stage esterification (depth <= 3)
            esterification_reactions = [
                "Esterification of Carboxylic Acids",
                "Transesterification",
                "Schotten-Baumann to ester",
                "O-alkylation of carboxylic acids with diazo compounds",
                "Acetic anhydride and alcohol to ester",
                "Oxidative esterification of primary alcohols"
            ]
            for r_name in esterification_reactions:
                if checker.check_reaction(r_name, rsmi):
                    if depth <= 3:
                        late_esterification = True
                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                        print(f"Found late-stage esterification at depth {depth}")
                    break

            # Check for Williamson ether synthesis or other ether formations
            ether_forming_reactions = [
                "Williamson Ether Synthesis",
                "{Williamson ether}",
                "Mitsunobu_phenole",
                "Chan-Lam etherification",
                "Alcohol to ether"
            ]
            ether_formation_found_in_step = False
            for r_name in ether_forming_reactions:
                if checker.check_reaction(r_name, rsmi):
                    ether_formations += 1
                    ether_formation_found_in_step = True
                    findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    print(f"Found ether formation at depth {depth}")
                    break

            # Check if this reaction involves aromatic fragments
            aromatic_reactants_in_ether_formation = 0
            if ether_formation_found_in_step:
                for reactant in reactants:
                    for ring in AROMATIC_RINGS_OF_INTEREST:
                        if checker.check_ring(ring, reactant):
                            aromatic_reactants_in_ether_formation += 1
                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                            print(f"Found aromatic reactant in ether formation: {reactant}")

                if aromatic_reactants_in_ether_formation >= 1:
                    has_aromatic_fragments += 1
                    print(f"Found aromatic fragment in ether formation at depth {depth}")

            # Check for polyether linker in product
            if checker.check_fg("Ether", product):
                ether_count_fg = len(checker.get_fg_atom_indices("Ether", product))
                if ether_count_fg >= 2:
                    has_polyether_linker = True
                    findings_json["atomic_checks"]["functional_groups"].append("Ether")
                    print(
                        f"Found polyether linker using FG check at depth {depth} with {ether_count_fg} ether groups"
                    )

            # Check for branching in synthesis (convergent indicator)
            if len(reactants) >= 2:
                complex_reactants = 0
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.GetNumAtoms() > 5:  # Non-trivial fragment
                        complex_reactants += 1

                if complex_reactants >= 2:
                    branching_points += 1
                    print(
                        f"Found branching point at depth {depth} with {complex_reactants} complex reactants"
                    )

        # Process children
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine if the strategy is present
    strategy_present = (
        (ether_formations >= 1 or has_polyether_linker)
        and (late_esterification or final_product_has_ester)
        and has_aromatic_fragments >= 1
        and branching_points >= 1
    )

    # Populate structural constraints based on flags
    if branching_points >= 1:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "convergent_step",
                "operator": ">=",
                "value": 1,
                "description": "At least one reaction step must be convergent, defined as having 2 or more reactants with more than 5 atoms."
            }
        })

    if has_aromatic_fragments >= 1:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "description": "At least one ether formation reaction must involve a specified aromatic ring in one of its reactants.",
                "min_occurrences": 1,
                "conditions": [
                    {
                        "target_type": "reaction",
                        "targets": [
                            "Williamson Ether Synthesis",
                            "{Williamson ether}",
                            "Mitsunobu_phenole",
                            "Chan-Lam etherification",
                            "Alcohol to ether"
                        ]
                    },
                    {
                        "target_type": "ring_system",
                        "targets": [
                            "benzene",
                            "naphthalene",
                            "anthracene",
                            "pyridine",
                            "furan",
                            "thiophene",
                            "pyrrole",
                            "indole",
                            "quinoline",
                            "isoquinoline"
                        ],
                        "scope": "reactant"
                    }
                ]
            }
        })

    if ether_formations >= 1 or has_polyether_linker:
        findings_json["structural_constraints"].append({
            "type": "alternative",
            "details": {
                "description": "The route must contain at least one ether formation reaction OR a product molecule with two or more ether groups.",
                "conditions": [
                    {
                        "type": "count",
                        "details": {
                            "target_type": "reaction",
                            "targets": [
                                "Williamson Ether Synthesis",
                                "{Williamson ether}",
                                "Mitsunobu_phenole",
                                "Chan-Lam etherification",
                                "Alcohol to ether"
                            ],
                            "operator": ">=",
                            "value": 1
                        }
                    },
                    {
                        "type": "count",
                        "details": {
                            "target_type": "functional_group",
                            "targets": [
                                "Ether"
                            ],
                            "scope": "product_molecule",
                            "operator": ">=",
                            "value": 2
                        }
                    }
                ]
            }
        })

    if late_esterification or final_product_has_ester:
        findings_json["structural_constraints"].append({
            "type": "alternative",
            "details": {
                "description": "The route must contain a late-stage esterification (depth <= 3) OR the final product must contain an ester group.",
                "conditions": [
                    {
                        "type": "positional",
                        "details": {
                            "target_type": "reaction",
                            "targets": [
                                "Esterification of Carboxylic Acids",
                                "Transesterification",
                                "Schotten-Baumann to ester",
                                "O-alkylation of carboxylic acids with diazo compounds",
                                "Acetic anhydride and alcohol to ester",
                                "Oxidative esterification of primary alcohols"
                            ],
                            "position": "depth_le_3"
                        }
                    },
                    {
                        "type": "positional",
                        "details": {
                            "target_type": "functional_group",
                            "targets": [
                                "Ester"
                            ],
                            "position": "final_product"
                        }
                    }
                ]
            }
        })

    print(f"Strategy detection results:")
    print(f"- Ether formations: {ether_formations}")
    print(f"- Late esterification: {late_esterification}")
    print(f"- Final product has ester: {final_product_has_ester}")
    print(f"- Polyether linker: {has_polyether_linker}")
    print(f"- Aromatic fragments: {has_aromatic_fragments}")
    print(f"- Branching points: {branching_points}")
    print(f"- Strategy present: {strategy_present}")

    return strategy_present, findings_json

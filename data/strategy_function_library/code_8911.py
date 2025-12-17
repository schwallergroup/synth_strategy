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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis with late-stage amide coupling between a carboxylic acid and an amine,
    where one fragment contains a piperidine ring and the other fragment undergoes nitrile reduction.
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

    # Initialize flags to track key features
    has_amide_coupling_final_step = False
    has_nitrile_reduction = False
    has_piperidine_fragment = False

    # Track which fragments have the required features
    piperidine_fragments = set()
    nitrile_reduction_fragments = set()

    # Track products of nitrile reduction
    nitrile_reduction_products = set()

    # Track final amide coupling reactants
    final_amide_reactants = []

    def dfs_traverse(node, depth=0, branch_id=None):
        nonlocal has_amide_coupling_final_step, has_nitrile_reduction, has_piperidine_fragment
        nonlocal piperidine_fragments, nitrile_reduction_fragments, nitrile_reduction_products, final_amide_reactants
        nonlocal findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for piperidine ring in molecule
            if checker.check_ring("piperidine", mol_smiles):
                has_piperidine_fragment = True
                findings_json["atomic_checks"]["ring_systems"].append("piperidine")
                if branch_id is not None:
                    piperidine_fragments.add(branch_id)

        elif node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for amide coupling in the final step
            if depth == 0:
                # Check for amide coupling reaction
                amide_coupling_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Carboxylic acid with primary amine to amide",
                    "Schotten-Baumann_amide"
                ]
                for reaction_name in amide_coupling_reactions:
                    if checker.check_reaction(reaction_name, rsmi):
                        # Verify one reactant has carboxylic acid and another has amine
                        acid_reactants = [
                            r for r in reactants_smiles if checker.check_fg("Carboxylic acid", r)
                        ]
                        amine_reactants = [
                            r
                            for r in reactants_smiles
                            if checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                        ]

                        if acid_reactants and amine_reactants and len(reactants_smiles) >= 2:
                            has_amide_coupling_final_step = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": [
                                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                                        "Carboxylic acid with primary amine to amide",
                                        "Schotten-Baumann_amide"
                                    ],
                                    "position": "last_stage"
                                }
                            })
                            if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                            if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"] and "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                                if any(checker.check_fg("Primary amine", r) for r in amine_reactants):
                                    findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                                elif any(checker.check_fg("Secondary amine", r) for r in amine_reactants):
                                    findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

                            final_amide_reactants = reactants_smiles
                            break # Found an amide coupling reaction

            # Check for nitrile reduction
            if checker.check_reaction("Reduction of nitrile to amine", rsmi):
                has_nitrile_reduction = True
                findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitrile to amine")
                if branch_id is not None:
                    nitrile_reduction_fragments.add(branch_id)
                # Track the product of nitrile reduction
                nitrile_reduction_products.add(product_smiles)

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # Only increase depth when going from chemical to reaction
            next_depth = depth + 1

        # Traverse children
        for i, child in enumerate(node.get("children", [])):
            # Create a unique branch ID for tracking fragments
            new_branch_id = f"{branch_id}_{i}" if branch_id is not None else str(i)
            dfs_traverse(child, next_depth, new_branch_id)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the final amide coupling connects a piperidine fragment with a product of nitrile reduction
    piperidine_in_final_reactants = False
    amine_from_nitrile_in_final_reactants = False

    if final_amide_reactants:
        piperidine_in_final_reactants = any(
            checker.check_ring("piperidine", r) for r in final_amide_reactants
        )
        if piperidine_in_final_reactants and "piperidine" not in findings_json["atomic_checks"]["ring_systems"]:
            findings_json["atomic_checks"]["ring_systems"].append("piperidine")

        amine_reactants = [
            r
            for r in final_amide_reactants
            if checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
        ]
        if amine_reactants:
            if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"] and any(checker.check_fg("Primary amine", r) for r in amine_reactants):
                findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
            if "Secondary amine" not in findings_json["atomic_checks"]["functional_groups"] and any(checker.check_fg("Secondary amine", r) for r in amine_reactants):
                findings_json["atomic_checks"]["functional_groups"].append("Secondary amine")

        for amine in amine_reactants:
            if amine in nitrile_reduction_products:
                amine_from_nitrile_in_final_reactants = True
                break

    # Check if all required features are present
    basic_requirements = has_amide_coupling_final_step and has_piperidine_fragment

    # For true convergent synthesis, verify that the final amide coupling connects the piperidine fragment with a product of nitrile reduction
    result = False
    if basic_requirements:
        # Check if piperidine and nitrile reduction are in different branches
        different_branches = not piperidine_fragments.intersection(nitrile_reduction_fragments)
        if different_branches:
            findings_json["structural_constraints"].append({
                "type": "disjoint_pathways",
                "details": {
                    "targets": [
                        "piperidine",
                        "Reduction of nitrile to amine"
                    ]
                }
            })

        # Check if the final amide coupling connects the required fragments
        convergent_structure = (
            piperidine_in_final_reactants and amine_from_nitrile_in_final_reactants
        )

        # The final check must be specific to the convergent strategy
        if convergent_structure and different_branches:
            result = True
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "reactant_with_carboxylic_acid",
                        "reactant_with_amine",
                        "reactant_with_piperidine"
                    ],
                    "context": "last_stage_reaction"
                }
            })
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "first_event": "Reduction of nitrile to amine",
                    "second_event": "last_stage_amide_coupling",
                    "relationship": "product_is_reactant"
                }
            })

    return result, findings_json
from typing import Tuple, Dict, List
import copy
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
    Detects a convergent synthesis strategy where one or more amide couplings are used to prepare a fragment,
    which is then used in a subsequent multi-component reaction. A convergent step is identified if it combines
    a fragment produced by a prior amide coupling, or if it combines at least two fragments that both contain amide groups.
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

    amide_coupling_count = 0
    amide_coupling_depths = []
    has_convergent_step = False
    amide_products = {}  # product SMILES -> depth

    # Define the full list of amide coupling reaction names for easy lookup
    AMIDE_COUPLING_REACTIONS = [
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Carboxylic acid with primary amine to amide",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Acyl chloride with secondary amine to amide",
        "Acylation of primary amines",
        "Acylation of secondary amines",
        "{Schotten-Baumann_amide}"
    ]

    def dfs_traverse(node, depth=0):
        nonlocal amide_coupling_count, amide_coupling_depths, has_convergent_step, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for amide coupling reactions using checker functions
            is_amide_coupling = False
            for reaction_name in AMIDE_COUPLING_REACTIONS:
                if checker.check_reaction(reaction_name, rsmi):
                    is_amide_coupling = True
                    if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                    break

            if is_amide_coupling:
                amide_coupling_count += 1
                amide_coupling_depths.append(depth)
                amide_products[product_smiles] = depth
                # print(f"Found amide coupling at depth {depth}, product: {product_smiles[:20]}...")

            # Check for convergent step (multiple reactants being combined)
            if len(reactants_smiles) >= 2:
                # Check if any reactant was previously formed by amide coupling
                amide_fragment_reactants = [r for r in reactants_smiles if r in amide_products]

                # Also check if any reactant contains an amide group
                amide_containing_reactants = []
                for r in reactants_smiles:
                    found_amide_fg = False
                    if checker.check_fg("Primary amide", r):
                        if "Primary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amide")
                        found_amide_fg = True
                    if checker.check_fg("Secondary amide", r):
                        if "Secondary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Secondary amide")
                        found_amide_fg = True
                    if checker.check_fg("Tertiary amide", r):
                        if "Tertiary amide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Tertiary amide")
                        found_amide_fg = True
                    if found_amide_fg:
                        amide_containing_reactants.append(r)

                if len(amide_containing_reactants) >= 1 and (
                    # Either a previously formed amide fragment is used
                    amide_fragment_reactants
                    or
                    # Or multiple amide-containing reactants are combined
                    len(amide_containing_reactants) >= 2
                ):
                    has_convergent_step = True
                    # print(f"Found convergent step using amide-containing fragment at depth {depth}")
                    # Add structural constraint for co-occurrence if multiple amide-containing reactants
                    if len(amide_containing_reactants) >= 2:
                        constraint_obj = {
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "amide_functional_group",
                                    "amide_functional_group"
                                ],
                                "scope": "reaction_reactants"
                            }
                        }
                        if constraint_obj not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(constraint_obj)

        # Process children (reactants in retrosynthesis)
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Determine if this matches our strategy
    strategy_detected = (
        amide_coupling_count >= 1  # At least one amide coupling
        and has_convergent_step  # And a convergent step
    )

    # Add structural constraints based on final flags
    if amide_coupling_count >= 1:
        constraint_obj = {
            "type": "count",
            "details": {
                "target": "amide_coupling_reaction",
                "operator": ">=",
                "value": 1
            }
        }
        if constraint_obj not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(constraint_obj)

    # The 'sequence' constraint is implicitly covered by 'has_convergent_step' being true
    # and the logic within dfs_traverse that checks for amide_fragment_reactants.
    # If has_convergent_step is true and it was due to amide_fragment_reactants, then the sequence is met.
    # We add it if both conditions are met.
    if has_convergent_step and amide_coupling_count >= 1:
        # This is a simplified representation, as the actual sequence check is more complex
        # and involves tracking which specific amide coupling product was used in a convergent step.
        # For this refactoring, we assume if both flags are true, the sequence is conceptually met.
        constraint_obj = {
            "type": "sequence",
            "details": {
                "before": "amide_coupling_reaction",
                "after": "multi_component_reaction_using_amide_product"
            }
        }
        if constraint_obj not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(constraint_obj)

    # print(f"Amide coupling strategy detection result: {strategy_detected}")
    # print(f"Amide coupling count: {amide_coupling_count}, at depths: {amide_coupling_depths}")
    # print(f"Has convergent step: {has_convergent_step}")

    return strategy_detected, findings_json

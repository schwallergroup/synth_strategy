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


CYCLIZATION_REACTIONS_OF_INTEREST = [
    "Formation of NOS Heterocycles",
    "Intramolecular transesterification/Lactone formation",
    "Paal-Knorr pyrrole synthesis",
    "Williamson Ether Synthesis (intra to epoxy)",
    "Production of 2H,1-benzopyrans",
    "Benzoxazole formation (intramolecular)",
    "Mitsunobu aryl ether (intramolecular)",
    "Intramolecular amination of azidobiphenyls (heterocycle formation)",
    "Intramolecular amination (heterocycle formation)",
    "Fischer indole",
    "Paal-Knorr pyrrole",
    "benzimidazole_derivatives_carboxylic-acid/ester",
    "benzimidazole_derivatives_aldehyde",
    "benzothiazole",
    "benzoxazole_arom-aldehyde",
    "benzoxazole_carboxylic-acid",
    "thiazole",
    "Pictet-Spengler",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a linear synthesis strategy that culminates in a late-stage
    cyclization step, identified by checking if the reaction matches any entry in the
    CYCLIZATION_REACTIONS_OF_INTEREST list.
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

    # Track reactions and their depths
    reactions_by_depth = {}

    # Track if we found a late cyclization
    found_late_cyclization = False

    # Track maximum depth to determine what "late-stage" means
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal reactions_by_depth, found_late_cyclization, max_depth, findings_json

        # Update max depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = [r for r in rsmi.split(">")[0].split(".") if r]
                product_smiles = rsmi.split(">")[-1]

                # Store reaction at this depth
                reactions_by_depth[depth] = {
                    "rsmi": rsmi,
                    "reactants": reactants_smiles,
                    "product": product_smiles,
                }

                ring_formed = False

                # Check for specific, named cyclization reactions
                for rxn in CYCLIZATION_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(rxn, rsmi):
                        ring_formed = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        break

                # Check if this is a late-stage step (depth 0 or 1)
                if ring_formed and depth <= 1:
                    found_late_cyclization = True
                    # Add the positional constraint to findings_json
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "ring_formation",
                            "position": "late_stage",
                            "description": "A cyclization reaction must occur in the last or second-to-last step of the synthesis (depth <= 1)."
                        }
                    })
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children with incremented depth based on node type
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have a linear synthesis (no convergent steps)
    is_linear = True
    for depth, reaction in reactions_by_depth.items():
        reactants = reaction["reactants"]
        # If any reaction has more than 2 non-empty reactants, it might be convergent
        if len([r for r in reactants if r]) > 2:
            is_linear = False
            print(f"Non-linear step detected at depth {depth} with {len(reactants)} reactants")
            # This condition means the negation constraint is NOT met, so we don't add it to findings_json
            break
    
    if is_linear:
        # If it is linear, add the negation constraint to findings_json
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "convergent_reaction_step",
                "description": "The route must be linear, meaning no reaction step has more than 2 reactants."
            }
        })

    # Determine if the cyclization is truly "late-stage" relative to the total synthesis
    late_stage_threshold = min(1, max_depth // 3)  # Consider the first third as "late-stage"

    print(f"Max depth: {max_depth}, Late-stage threshold: {late_stage_threshold}")
    print(f"Is linear: {is_linear}, Found late cyclization: {found_late_cyclization}")

    # Return True if we have a linear synthesis with late cyclization
    result = is_linear and found_late_cyclization
    return result, findings_json

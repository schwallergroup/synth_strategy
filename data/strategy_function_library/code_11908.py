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


# List of heteroaromatic rings to check
HETEROAROMATIC_RINGS_OF_INTEREST = [
    "pyridine", "pyrazole", "imidazole", "oxazole", "thiazole",
    "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole",
    "indole", "quinoline", "isoquinoline", "purine", "benzoxazole",
    "benzothiazole", "benzimidazole", "furan", "thiophene", "isoxazole",
    "isothiazole", "oxadiazole", "thiadiazole",
]

# List of common coupling reactions
COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki", "Buchwald-Hartwig", "Stille", "Negishi", "Heck",
    "Sonogashira", "Ullmann-Goldberg", "N-arylation", "Chan-Lam",
    "Kumada", "Hiyama-Denmark", "Catellani", "decarboxylative_coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage functionalization strategy where a heteroaromatic group is introduced
    in the final synthetic step via a specified cross-coupling reaction. This is confirmed by
    checking that: (1) the reaction occurs at the final step (depth=1), (2) it is one of the
    enumerated coupling reactions in COUPLING_REACTIONS_OF_INTEREST, and (3) a heteroaromatic
    ring from HETEROAROMATIC_RINGS_OF_INTEREST is present in one of the reactants and also
    in the product.
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

    final_step_heteroaryl_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_heteroaryl_coupling, findings_json

        # Skip leaf nodes (starting materials)
        if node["type"] == "mol" and node.get("in_stock", False):
            return

        if node["type"] == "reaction" and depth == 1:  # Final step only
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains a heteroaromatic ring
                product_has_heteroaromatic = False
                product_heteroaromatic_rings = []

                for ring in HETEROAROMATIC_RINGS_OF_INTEREST:
                    if checker.check_ring(ring, product):
                        product_has_heteroaromatic = True
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                        product_heteroaromatic_rings.append(ring)

                if not product_has_heteroaromatic:
                    return

                # Check which reactants contain heteroaromatic rings
                reactant_with_heteroaromatic = None
                reactant_heteroaromatic_rings = []

                for i, reactant in enumerate(reactants):
                    reactant_rings = []
                    for ring in HETEROAROMATIC_RINGS_OF_INTEREST:
                        if checker.check_ring(ring, reactant):
                            reactant_rings.append(ring)
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)

                    if reactant_rings:
                        reactant_with_heteroaromatic = i
                        reactant_heteroaromatic_rings = reactant_rings

                # Check if this is a coupling reaction
                is_coupling_reaction = False
                detected_coupling_reaction = None
                for reaction_type in COUPLING_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_coupling_reaction = True
                        detected_coupling_reaction = reaction_type
                        if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

                # Determine if this is a late-stage heteroaryl introduction
                if (
                    product_has_heteroaromatic
                    and reactant_with_heteroaromatic is not None
                    and is_coupling_reaction
                ):
                    final_step_heteroaryl_coupling = True
                    # Add the structural constraint if all conditions are met
                    findings_json["structural_constraints"].append({
                        "type": "co-occurrence",
                        "details": {
                            "targets": [
                                "any_from_COUPLING_REACTIONS_OF_INTEREST",
                                "any_from_HETEROAROMATIC_RINGS_OF_INTEREST_in_product",
                                "any_from_HETEROAROMATIC_RINGS_OF_INTEREST_in_reactant"
                            ],
                            "scope": "single_reaction",
                            "positional_constraint": "last_stage"
                        }
                    })

            except Exception as e:
                pass

        # Traverse children with modified depth calculation
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return final_step_heteroaryl_coupling, findings_json

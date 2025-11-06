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


CC_COUPLING_REACTIONS = [
    "Suzuki",
    "Negishi",
    "Stille",
    "Kumada",
    "Heck",
    "Sonogashira",
    "Catellani",
    "Hiyama-Denmark",
    "Aryllithium cross-coupling",
    "decarboxylative_coupling",
    "Grignard",
    "Wittig",
    "Aldol condensation",
    "Knoevenagel Condensation",
    "Michael addition",
    "Diels-Alder",
    "Friedel-Crafts alkylation",
    "Friedel-Crafts acylation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects linear synthetic routes of at least two steps that feature a late-stage carbon-carbon bond formation involving at least one aromatic reactant. The reaction is identified if its type is present in the `CC_COUPLING_REACTIONS` list.
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

    # Track if we found a late-stage C-C bond formation
    late_stage_cc_bond = False
    # Track if the synthesis is linear (no branching)
    is_linear = True
    # Track reaction depths
    max_depth = 0
    # List to store all reaction nodes with their depths
    reaction_nodes = []

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_cc_bond, is_linear, max_depth, findings_json

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "reaction":
            # Store reaction node with its depth for later analysis
            reaction_nodes.append((node, depth))

            # Check if this is a branching point (convergent synthesis)
            mol_children = [child for child in node.get("children", []) if child["type"] == "mol"]
            if len(mol_children) > 2:
                is_linear = False
                # Record the structural constraint for negation of convergent step
                findings_json["structural_constraints"].append({
                    "type": "negation",
                    "details": {
                        "target": "convergent_step",
                        "definition": "A reaction with more than 2 molecule reactants."
                    }
                })

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Define what "late stage" means based on max_depth
    late_stage_threshold = max(1, max_depth // 3)

    # Check reactions for C-C bond formation between aromatic fragments
    for node, depth in reaction_nodes:
        # Only check reactions in the late stage
        if depth <= late_stage_threshold:
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]

                # Check if this is a known C-C coupling reaction
                cc_coupling_reaction = False
                for rxn_type in CC_COUPLING_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        cc_coupling_reaction = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                if cc_coupling_reaction:
                    reactants = reactants_part.split(".")

                    # Check if reactants contain aromatic rings
                    aromatic_reactants_found = False
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                for atom in mol.GetAtoms():
                                    if atom.GetIsAromatic():
                                        aromatic_reactants_found = True
                                        findings_json["atomic_checks"]["ring_systems"].append("aromatic_ring")
                                        break
                            if aromatic_reactants_found:
                                break
                        except Exception:
                            continue

                    if aromatic_reactants_found:
                        late_stage_cc_bond = True
                        # Record the co-occurrence structural constraint
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "scope": "same_reaction_step",
                                "targets": [
                                    "cc_coupling_reaction",
                                    "aromatic_reactant"
                                ]
                            }
                        })
                        # Record the positional structural constraint
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "cc_coupling_reaction",
                                "position": "late_stage_1/3"
                            }
                        })

    # A linear synthesis with late-stage C-C bond formation
    result = is_linear and late_stage_cc_bond and max_depth >= 3

    # Record the count structural constraint if met
    if max_depth >= 3:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction_steps",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json

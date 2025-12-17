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


MULTI_COMPONENT_COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl halide",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Heck terminal vinyl",
    "Negishi coupling",
    "Stille reaction_aryl",
    "Ugi reaction",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "A3 coupling",
    "Chan-Lam amine",
    "Chan-Lam alcohol",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage (final 2 steps), multi-component (>= 3 reactants) reaction. The reaction must be one of the types specified in the MULTI_COMPONENT_COUPLING_REACTIONS list.
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

    has_multi_component_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_multi_component_coupling, findings_json

        if has_multi_component_coupling:
            return

        if node["type"] == "reaction" and depth <= 2:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                valid_reactants = [r for r in reactants_smiles if r.strip()]

                if len(valid_reactants) >= 3:
                    # Record the reactant count constraint if met
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "reactants",
                            "operator": ">=",
                            "value": 3
                        }
                    })
                    # Check if it's one of the specified multi-component coupling reactions
                    for reaction_type in MULTI_COMPONENT_COUPLING_REACTIONS:
                        if checker.check_reaction(reaction_type, rsmi):
                            has_multi_component_coupling = True
                            # Record the named reaction finding
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                            # Record the positional constraint if met
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": "multi-component coupling reaction",
                                    "position": "within_last_2_stages"
                                }
                            })
                            break
            except Exception:
                # Silently ignore malformed reaction nodes
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            # Depth remains the same when traversing from reaction to chemical
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return has_multi_component_coupling, findings_json

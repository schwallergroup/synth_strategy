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


FURAN_LIKE_HETEROCYCLES = ["furan", "oxazole", "isoxazole", "thiazole", "oxadiazole"]
COUPLING_REACTIONS = [
    "Suzuki",
    "Negishi",
    "Stille",
    "Heck",
    "Sonogashira",
    "Buchwald-Hartwig",
    "Ullmann",
    "N-arylation",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Aryllithium cross-coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage (final 3 steps) coupling reaction that joins a pyridine-containing fragment with a fragment containing a furan-like heterocycle. This is identified by checking for specific named coupling reactions where the product contains both heterocyclic systems, which are also present across the reactant pool. The specific heterocycles and coupling reactions are defined in the module-level constants `FURAN_LIKE_HETEROCYCLES` and `COUPLING_REACTIONS`.
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

    found_heterocycle_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycle_coupling, findings_json

        if found_heterocycle_coupling:
            return

        # For reaction nodes, check for coupling reactions
        if node["type"] == "reaction" and depth <= 2:  # Late stage (depth 0, 1, or 2)
            # Record positional constraint if met
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "any_coupling_reaction",
                    "position": "depth <= 2"
                }
            })

            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a coupling reaction
                is_coupling = False
                for rxn_type in COUPLING_REACTIONS:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_coupling = True
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)

                # Check for C-C or C-N bond formation between heterocycles
                if is_coupling:
                    # Check for heterocycles in reactants
                    furan_like_in_reactants = False
                    reactant_heterocycles_found = []
                    for het in FURAN_LIKE_HETEROCYCLES:
                        if any(checker.check_ring(het, r) for r in reactants):
                            furan_like_in_reactants = True
                            if het not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(het)
                            reactant_heterocycles_found.append(het)

                    pyridine_in_reactants = False
                    if any(
                        checker.check_ring("pyridine", r) for r in reactants
                    ):
                        pyridine_in_reactants = True
                        if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("pyridine")
                        reactant_heterocycles_found.append("pyridine")

                    # Check if product has both heterocycles
                    product_has_furan_like = False
                    product_heterocycles_found = []
                    for het in FURAN_LIKE_HETEROCYCLES:
                        if checker.check_ring(het, product):
                            product_has_furan_like = True
                            if het not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(het)
                            product_heterocycles_found.append(het)

                    product_has_pyridine = False
                    if checker.check_ring("pyridine", product):
                        product_has_pyridine = True
                        if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("pyridine")
                        product_heterocycles_found.append("pyridine")

                    # If product has both heterocycles and at least one was in reactants,
                    # this is likely a heterocycle coupling
                    if (
                        product_has_furan_like
                        and product_has_pyridine
                    ):
                        # Record co-occurrence constraint if met
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "scope": "reaction_product",
                                "targets": [
                                    "pyridine",
                                    "any_furan_like_heterocycle"
                                ]
                            }
                        })

                        if (furan_like_in_reactants or pyridine_in_reactants):
                            found_heterocycle_coupling = True
                            # Record count constraint if met
                            findings_json["structural_constraints"].append({
                                "type": "count",
                                "details": {
                                    "scope": "reaction_reactants",
                                    "target_set": [
                                        "pyridine",
                                        "any_furan_like_heterocycle"
                                    ],
                                    "operator": ">=",
                                    "value": 1
                                }
                            })

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                # From reaction to chemical, depth remains the same
                dfs_traverse(child, depth)
            else:
                # From chemical to reaction, depth increases
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_heterocycle_coupling, findings_json

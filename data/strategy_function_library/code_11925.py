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


HETEROCYCLES_OF_INTEREST = ["pyridine", "pyrimidine", "pyrazine", "pyridazine", "quinoline", "isoquinoline"]
FUNCTIONAL_GROUPS_FOR_MODIFICATION = [
    "Primary alcohol", "Secondary alcohol", "Tertiary alcohol", "Ketone", "Aldehyde",
    "Primary amine", "Secondary amine", "Tertiary amine", "Primary amide", "Secondary amide",
    "Tertiary amide", "Primary halide", "Secondary halide", "Tertiary halide",
    "Aromatic halide", "Nitrile", "Carboxylic acid", "Ester", "Nitro group",
    "Sulfonamide", "Sulfone", "Phosphate ester",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a linear synthesis strategy that maintains a specific heterocyclic scaffold
    throughout the route, while performing sequential functional group modifications.
    The scaffolds of interest are defined in the HETEROCYCLES_OF_INTEREST list, and the
    modifications are tracked by checking for changes in the functional groups defined
    in the FUNCTIONAL_GROUPS_FOR_MODIFICATION list. The strategy is flagged if a
    scaffold is seen and maintained, the synthesis is linear (i.e., no convergent
    steps involving two scaffold-containing molecules), and at least two modification
    steps occur.
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

    scaffold_maintained = True
    fg_modifications = 0
    is_linear = True
    seen_heterocycle = False

    # Define the structural constraints from the JSON metadata
    structural_constraints_metadata = [
        {
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "target_scaffold_presence"
                ],
                "description": "At least one of the specified heterocyclic scaffolds must be present somewhere in the synthesis route."
            }
        },
        {
            "type": "negation",
            "details": {
                "target": "ring_destruction",
                "description": "A target heterocyclic scaffold, once present, must not be destroyed in any subsequent reaction."
            }
        },
        {
            "type": "count",
            "details": {
                "target": "scaffold_containing_reactants",
                "operator": "<=",
                "value": 1,
                "scope": "per_reaction",
                "description": "Ensures linearity by allowing at most one scaffold-containing reactant per reaction step."
            }
        },
        {
            "type": "count",
            "details": {
                "target": "functional_group_modification_on_scaffold",
                "operator": ">=",
                "value": 2,
                "description": "At least two modifications involving the specified functional groups must occur on the scaffold-containing molecule."
            }
        }
    ]

    def dfs_traverse(node, depth=0):
        nonlocal scaffold_maintained, fg_modifications, is_linear, seen_heterocycle, findings_json

        new_depth = depth
        if node["type"] != "reaction":
            new_depth = depth + 1

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            main_reactants = []
            for r in reactants:
                if not r:
                    continue
                mol = Chem.MolFromSmiles(r)
                if mol and mol.GetNumHeavyAtoms() >= 5:
                    main_reactants.append(r)

            reactant_has_heterocycle = False
            heterocycle_reactants = []

            for r in main_reactants:
                for ring in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(ring, r):
                        reactant_has_heterocycle = True
                        heterocycle_reactants.append(r)
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                        seen_heterocycle = True
                        break

            product_has_heterocycle = False
            if product:
                for ring in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(ring, product):
                        product_has_heterocycle = True
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                        seen_heterocycle = True
                        break

            if len(heterocycle_reactants) > 1:
                is_linear = False
                # Record the violation of linearity constraint
                for constraint in structural_constraints_metadata:
                    if constraint["type"] == "count" and constraint["details"]["target"] == "scaffold_containing_reactants":
                        if constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(constraint)

            if reactant_has_heterocycle and not product_has_heterocycle:
                scaffold_maintained = False
                # Record the ring destruction
                if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")
                for constraint in structural_constraints_metadata:
                    if constraint["type"] == "negation" and constraint["details"]["target"] == "ring_destruction":
                        if constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(constraint)

            if reactant_has_heterocycle and product_has_heterocycle:
                for reactant in heterocycle_reactants:
                    for fg in FUNCTIONAL_GROUPS_FOR_MODIFICATION:
                        reactant_has_fg = checker.check_fg(fg, reactant)
                        product_has_fg = checker.check_fg(fg, product)

                        if (reactant_has_fg and not product_has_fg) or \
                           (not reactant_has_fg and product_has_fg):
                            fg_modifications += 1
                            if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(fg)
                            break

        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    strategy_present = (
        seen_heterocycle and scaffold_maintained and is_linear and fg_modifications >= 2
    )

    # Record structural constraints based on final flags
    if seen_heterocycle:
        for constraint in structural_constraints_metadata:
            if constraint["type"] == "co-occurrence" and "target_scaffold_presence" in constraint["details"]["targets"]:
                if constraint not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(constraint)

    if fg_modifications >= 2:
        for constraint in structural_constraints_metadata:
            if constraint["type"] == "count" and constraint["details"]["target"] == "functional_group_modification_on_scaffold":
                if constraint not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(constraint)

    return strategy_present, findings_json

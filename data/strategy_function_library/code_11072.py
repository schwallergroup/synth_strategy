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


# List of heterocycles to check for incorporation
HETEROCYCLES_FOR_INCORPORATION = [
    "imidazole", "pyrrole", "pyridine", "pyrazole", "oxazole",
    "thiazole", "furan", "thiophene", "pyrimidine", "pyrazine",
    "indole", "benzimidazole",
]

# Coupling reactions that are strong indicators of convergent fragment coupling
COUPLING_REACTIONS_FOR_INCORPORATION = [
    "Suzuki", "Negishi", "Stille", "Heck", "Buchwald-Hartwig",
    "N-arylation", "Sonogashira", "Ullmann-Goldberg",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the late-stage (final two steps), convergent incorporation of a pre-formed heterocycle via a specified coupling reaction. This strategy is identified if a reaction is one of the types listed in `COUPLING_REACTIONS_FOR_INCORPORATION`, one reactant contains a heterocycle from `HETEROCYCLES_FOR_INCORPORATION`, this heterocycle is present in the product, and at least one other reactant is a complex fragment (>= 8 atoms).
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

    has_heterocycle_incorporation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_heterocycle_incorporation, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for heterocycle incorporation in late-stage reactions (depth <= 2)
            if depth <= 2:
                # Add positional constraint if met
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "convergent_heterocycle_incorporation",
                        "position": "within_final_2_steps"
                    }
                })

                # Check if this is a coupling reaction
                is_coupling = False
                detected_coupling_reaction = None
                for reaction_type in COUPLING_REACTIONS_FOR_INCORPORATION:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_coupling = True
                        detected_coupling_reaction = reaction_type
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

                if is_coupling:
                    # Check if one reactant has a heterocycle
                    heterocycle_reactant_idx = None
                    heterocycle_found = None
                    reactant_contains_heterocycle = False

                    for idx, reactant in enumerate(reactants):
                        for heterocycle in HETEROCYCLES_FOR_INCORPORATION:
                            if checker.check_ring(heterocycle, reactant):
                                heterocycle_reactant_idx = idx
                                heterocycle_found = heterocycle
                                reactant_contains_heterocycle = True
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                                break
                        if heterocycle_reactant_idx is not None:
                            break

                    # Check if product has the same heterocycle
                    product_contains_same_heterocycle = False
                    if heterocycle_reactant_idx is not None and checker.check_ring(
                        heterocycle_found, product
                    ):
                        product_contains_same_heterocycle = True
                        # The heterocycle is already added to ring_systems when found in reactant

                    if is_coupling and reactant_contains_heterocycle and product_contains_same_heterocycle:
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "named_coupling_reaction",
                                    "reactant_contains_heterocycle",
                                    "product_contains_same_heterocycle"
                                ]
                            }
                        })

                        # Verify the other reactant is complex enough (convergent synthesis)
                        other_reactants = [
                            r for i, r in enumerate(reactants) if i != heterocycle_reactant_idx
                        ]

                        for other_reactant in other_reactants:
                            other_mol = Chem.MolFromSmiles(other_reactant)
                            if other_mol and other_mol.GetNumAtoms() >= 8:  # Complex enough
                                findings_json["structural_constraints"].append({
                                    "type": "count",
                                    "details": {
                                        "target": "atoms_in_coupling_partner",
                                        "operator": ">=",
                                        "value": 8
                                    }
                                })
                                has_heterocycle_incorporation = True
                                return  # Found what we're looking for

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return has_heterocycle_incorporation, findings_json

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


# Module-level constants for enumeration
NITROGEN_HETEROCYCLES = [
    "pyrazine", "pyrazole", "piperidine", "piperazine", "pyridine", "imidazole",
    "triazole", "tetrazole", "morpholine", "thiomorpholine", "indole", "quinoline",
    "isoquinoline", "purine", "carbazole", "acridine", "benzimidazole",
    "benzotriazole", "benzoxazole", "benzothiazole", "pyrrolidine", "aziridine",
    "azetidine", "azepane",
]

COUPLING_REACTIONS = [
    "Suzuki", "Buchwald-Hartwig", "N-arylation", "Stille", "Negishi", "Heck",
    "Sonogashira", "Ullmann-Goldberg", "Kumada", "Hiyama-Denmark", "N-arylation_heterocycles",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthesis strategy involving coupling of multiple nitrogen heterocycles.
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

    # Track if we found the pattern
    found_n_heterocycle_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_n_heterocycle_coupling, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a coupling reaction
                is_coupling = False
                for reaction in COUPLING_REACTIONS:
                    if checker.check_reaction(reaction, rsmi):
                        is_coupling = True
                        if reaction not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction)
                        break

                if is_coupling:
                    # Count reactants that contain a heterocycle of interest
                    reactants_with_heterocycles = 0
                    found_heterocycles_in_reactants = []
                    for r_smiles in reactants_smiles:
                        for heterocycle in NITROGEN_HETEROCYCLES:
                            if checker.check_ring(heterocycle, r_smiles):
                                reactants_with_heterocycles += 1
                                if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                                found_heterocycles_in_reactants.append(heterocycle)
                                break  # Move to next reactant

                    # If at least two reactants have heterocycles, check the product
                    if reactants_with_heterocycles >= 2:
                        product_has_heterocycle = False
                        found_heterocycle_in_product = None
                        for heterocycle in NITROGEN_HETEROCYCLES:
                            if checker.check_ring(heterocycle, product_smiles):
                                product_has_heterocycle = True
                                if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                                found_heterocycle_in_product = heterocycle
                                break
                        
                        if product_has_heterocycle:
                            found_n_heterocycle_coupling = True
                            # Add structural constraint if all conditions are met
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "A coupling reaction from the specified list occurs.",
                                        "At least two reactants in the coupling reaction contain a specified nitrogen heterocycle.",
                                        "The product of the coupling reaction contains a specified nitrogen heterocycle."
                                    ],
                                    "scope": "single_reaction_step"
                                }
                            })

            except Exception:
                # Silently ignore errors in reaction processing
                pass

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_n_heterocycle_coupling, findings_json

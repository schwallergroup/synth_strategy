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


HETEROCYCLES_OF_INTEREST = [
    "indole", "benzimidazole", "benzoxazole", "benzothiazole", "pyrrole",
    "imidazole", "oxazole", "thiazole", "triazole", "tetrazole", "pyridine",
    "pyrimidine", "pyrazine", "pyridazine", "furan", "pyran", "dioxane",
    "tetrahydrofuran", "tetrahydropyran", "oxirane", "oxetane", "oxolane",
    "oxane", "dioxolane", "pyrrolidine", "piperidine", "piperazine",
    "morpholine", "thiomorpholine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a one-pot reaction where an azide-containing reactant undergoes cyclization to form a new heterocycle.
    This is confirmed by checking for an azide functional group in any reactant and the formation of a new ring in the product,
    where the ring structure is one of those specified in the HETEROCYCLES_OF_INTEREST list.
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

    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern, findings_json
        if found_pattern:  # Optimization: stop traversing if pattern is found
            return

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # This strategy requires an azide in at least one reactant.
                has_azide_reactant = False
                for r in reactants_smiles:
                    if checker.check_fg("Azide", r):
                        has_azide_reactant = True
                        if "Azide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Azide")
                        break

                if has_azide_reactant:
                    # Check for the formation of a new heterocycle from the predefined list.
                    for ring in HETEROCYCLES_OF_INTEREST:
                        is_in_product = checker.check_ring(ring, product_smiles)
                        is_in_reactants = any(checker.check_ring(ring, r) for r in reactants_smiles)
                        
                        if is_in_product and not is_in_reactants:
                            found_pattern = True
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)
                            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                            
                            # Add structural constraint if both conditions are met
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "Azide",
                                        "ring_formation"
                                    ],
                                    "scope": "reaction_step",
                                    "description": "An azide functional group must be present in a reactant that leads to the formation of a new heterocyclic ring in the same reaction step."
                                }
                            })
                            return # Pattern found, no need to check other rings or traverse children

            except Exception:
                # Silently ignore errors in reaction processing to ensure robustness.
                pass

        # Traverse children if the pattern has not yet been found.
        if not found_pattern:
            for child in node.get("children", []):
                # New logic for depth calculation:
                # Depth increases only when traversing from a chemical node to a reaction node.
                # Depth remains the same when traversing from a reaction node to a chemical node.
                if node["type"] == "reaction":
                    dfs_traverse(child, depth)
                else: # Assuming 'chemical' or other types
                    dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_pattern, findings_json
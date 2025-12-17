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
    "furan", "pyran", "dioxane", "tetrahydrofuran", "tetrahydropyran",
    "oxirane", "oxetane", "oxolane", "oxane", "dioxolane", "dioxolene",
    "pyrrole", "pyridine", "pyrazole", "imidazole", "oxazole", "thiazole",
    "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole",
    "pyrrolidine", "piperidine", "piperazine", "morpholine", "thiomorpholine",
    "aziridine", "azetidine", "thiophene", "thiopyran", "thiirane",
    "thietane", "isoxazole", "isothiazole", "oxadiazole", "thiadiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if heterocycle formation occurs early in the synthesis (high depth).
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

    heterocycle_formation_found = False
    formation_depth = -1
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_found, formation_depth, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            try:
                # Get product and reactants
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Check if a heterocycle is formed in this reaction
                reactant_smiles = reactants_part.split(".")

                # Check if product contains a heterocycle that none of the reactants contain
                product_heterocycles = set()
                for ring in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(ring, product_part):
                        product_heterocycles.add(ring)

                # Check which heterocycles are in reactants
                reactant_heterocycles = set()
                for reactant in reactant_smiles:
                    for ring in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(ring, reactant):
                            reactant_heterocycles.add(ring)

                # If there's a heterocycle in the product that's not in any reactant, it was formed
                new_heterocycles = product_heterocycles - reactant_heterocycles
                if new_heterocycles:
                    heterocycle_formation_found = True
                    formation_depth = depth
                    # Record the specific heterocycles found
                    for h_ring in new_heterocycles:
                        if h_ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(h_ring)
                    # Record that a ring formation occurred
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

            except Exception:
                # Silently ignore errors in reaction processing
                pass

        # Continue traversal
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemicals)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reactions)
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it "early" if it happens in the first half of the synthesis
    # In retrosynthetic analysis, high depth = early stage
    is_early = formation_depth >= max_depth / 2 if max_depth > 0 and formation_depth >= 0 else False

    result = heterocycle_formation_found and is_early

    if result:
        # Add the structural constraint if the overall condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "ring_formation",
                "position": "first_half"
            }
        })

    return result, findings_json

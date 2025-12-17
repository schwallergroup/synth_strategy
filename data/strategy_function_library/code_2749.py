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


# Define pairs of saturated and aromatic heterocycles to check
HETEROCYCLE_AROMATIZATION_PAIRS = [
    # Nitrogen-containing
    ("piperidine", "pyridine"),
    ("pyrrolidine", "pyrrole"),
    ("piperazine", "pyrazine"),
    # Oxygen-containing
    ("tetrahydrofuran", "furan"),
    ("tetrahydropyran", "pyran"),
    ("dioxane", "dioxene"),
    # Sulfur-containing
    ("thiomorpholine", "thiopyran"),
    ("thiolane", "thiophene"),
    ("thiane", "thiopyran"),
    # Mixed heterocycles
    ("oxazolidine", "oxazole"),
    ("thiazolidine", "thiazole"),
    ("imidazolidine", "imidazole"),
    ("pyrroline", "pyrrole"),
    ("dihydropyridine", "pyridine"),
    ("dihydropyrazine", "pyrazine"),
    ("dihydropyrimidine", "pyrimidine"),
]

# List of relevant dehydrogenation reaction types (Arene hydrogenation works in reverse)
DEHYDROGENATION_REACTION_TYPES = [
    "Arene hydrogenation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy involving late-stage heterocycle modification,
    specifically the conversion of a saturated heterocycle to an aromatic one.
    In retrosynthetic analysis, this means finding an aromatic heterocycle
    that is converted to a saturated one.
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

    found_late_stage_heterocycle_mod = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_heterocycle_mod, findings_json

        if node["type"] == "reaction" and depth <= 2:  # Check reactions up to depth 2 (late stage)
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # In retrosynthesis, the product of the reaction is the reactant in forward synthesis
                # and the reactants are the products in forward synthesis
                for reactant in reactants:
                    for sat_ring, arom_ring in HETEROCYCLE_AROMATIZATION_PAIRS:
                        # Check if fwd-product contains aromatic ring and fwd-reactant contains saturated ring
                        if checker.check_ring(arom_ring, reactant):
                            if arom_ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(arom_ring)
                        if checker.check_ring(sat_ring, product):
                            if sat_ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(sat_ring)

                        if checker.check_ring(arom_ring, reactant) and checker.check_ring(
                            sat_ring, product
                        ):
                            # Check for specific dehydrogenation reaction types
                            for rxn_type in DEHYDROGENATION_REACTION_TYPES:
                                if checker.check_reaction(rxn_type, rsmi):
                                    found_late_stage_heterocycle_mod = True
                                    if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                                    
                                    # Add structural constraints if found
                                    findings_json["structural_constraints"].append({
                                        "type": "positional",
                                        "details": {
                                            "target": "Arene hydrogenation",
                                            "position": "late_stage (depth <= 2)"
                                        }
                                    })
                                    findings_json["structural_constraints"].append({
                                        "type": "co-occurrence",
                                        "details": {
                                            "targets": [
                                                "Arene hydrogenation",
                                                "aromatic_ring_destruction",
                                                "saturated_ring_formation"
                                            ],
                                            "scope": "same_step"
                                        }
                                    })
                                    # Once found, no need to continue checking this reaction step
                                    return

            except Exception:
                # Silently ignore errors in reaction processing
                pass

        # Traverse children
        if not found_late_stage_heterocycle_mod:
            for child in node.get("children", []):
                # New logic for depth calculation
                new_depth = depth
                if node["type"] != "reaction": # If current node is chemical, depth increases
                    new_depth = depth + 1
                # If current node is reaction, depth remains the same
                dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return found_late_stage_heterocycle_mod, findings_json

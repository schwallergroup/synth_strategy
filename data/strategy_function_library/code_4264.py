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
    "furan", "pyran", "dioxane", "tetrahydrofuran", "tetrahydropyran", "oxirane",
    "oxetane", "oxolane", "oxane", "dioxolane", "dioxolene", "trioxane",
    "dioxepane", "pyrrole", "pyridine", "pyrazole", "imidazole", "oxazole",
    "thiazole", "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole",
    "pyrrolidine", "piperidine", "piperazine", "morpholine", "thiomorpholine",
    "aziridine", "azetidine", "azepane", "diazepane", "indole", "quinoline",
    "isoquinoline", "purine", "carbazole", "acridine", "thiophene", "thiopyran",
    "thiirane", "thietane", "thiolane", "thiane", "dithiane", "dithiolane",
    "benzothiophene", "oxathiolane", "dioxathiolane", "thiazolidine",
    "oxazolidine", "isoxazole", "isothiazole", "oxadiazole", "thiadiazole",
    "benzoxazole", "benzothiazole", "benzimidazole", "pteridin", "phenothiazine",
    "phenoxazine", "dibenzofuran", "dibenzothiophene", "xanthene",
    "thioxanthene", "pyrroline", "pyrrolidone", "imidazolidine", "porphyrin",
    "indazole", "benzotriazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves a late-stage formation (in the final 3 steps)
    of a specific heterocyclic ring. The list of heterocycles checked is defined in the
    module-level constant HETEROCYCLES_OF_INTEREST.
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

    ring_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ring_formation_detected, findings_json

        if node["type"] == "reaction" and depth <= 3:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Direct check for heterocycle formation
                for ring_name in HETEROCYCLES_OF_INTEREST:
                    # Check if product contains the heterocyclic ring
                    if checker.check_ring(ring_name, product_smiles):
                        # Check if any reactant already has this ring
                        reactant_has_ring = False
                        for r in reactants_smiles:
                            if checker.check_ring(ring_name, r):
                                reactant_has_ring = True
                                break

                        if not reactant_has_ring:
                            ring_formation_detected = True
                            if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                            # Add structural constraint if not already present
                            positional_constraint = {
                                "type": "positional",
                                "details": {
                                    "target": "ring_formation",
                                    "position": "within_final_4_steps"
                                }
                            }
                            if positional_constraint not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(positional_constraint)
                            return
            except Exception:
                # Silently ignore errors in reaction processing
                pass

        # Traverse children
        for child in node.get("children", []):
            if ring_formation_detected:
                return
            
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "chemical"
                dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return ring_formation_detected, findings_json
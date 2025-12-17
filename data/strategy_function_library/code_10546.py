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
    "pyridine", "pyrrole", "furan", "thiophene", "imidazole", "oxazole",
    "thiazole", "pyrazole", "isoxazole", "isothiazole", "triazole",
    "tetrazole", "pyrimidine", "pyrazine", "pyridazine", "indole",
    "benzofuran", "benzothiophene", "benzimidazole", "benzoxazole",
    "benzothiazole", "quinoline", "isoquinoline", "piperidine",
    "piperazine", "morpholine", "thiomorpholine",
]

COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Heck terminal vinyl",
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl halide",
    "Stille reaction_aryl",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Negishi coupling",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthesis route contains at least two instances of a specific coupling reaction that involves one of a predefined list of heterocycles. A reaction is counted if: 1. It matches a reaction type from the `COUPLING_REACTIONS_OF_INTEREST` list (e.g., Suzuki, Buchwald-Hartwig). 2. The product contains a heterocycle from the `HETEROCYCLES_OF_INTEREST` list that was not present in all of the reactants.
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

    heterocycle_couplings = 0

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_couplings, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    return

                parts = rsmi.split(">")
                if len(parts) < 3:
                    return

                reactants = parts[0].split(".")
                product = parts[-1]

                is_coupling = False
                for reaction_type in COUPLING_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_coupling = True
                        if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

                if is_coupling:
                    product_heterocycles = []
                    for ring in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(ring, product):
                            product_heterocycles.append(ring)
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)

                    new_heterocycle = False
                    for ring in product_heterocycles:
                        if not all(checker.check_ring(ring, r) for r in reactants):
                            new_heterocycle = True
                            break

                    if new_heterocycle:
                        heterocycle_couplings += 1
            except Exception:
                # Silently ignore errors in single reaction processing
                pass

        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = heterocycle_couplings >= 2
    if result:
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "coupling_reaction_introducing_a_heterocycle",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json

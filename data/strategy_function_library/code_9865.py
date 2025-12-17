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
    "furan",
    "pyrrole",
    "thiophene",
    "pyrazole",
    "imidazole",
    "oxazole",
    "thiazole",
    "isoxazole",
    "isothiazole",
    "triazole",
    "tetrazole",
    "pyridine",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
]

COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki",
    "Negishi",
    "Stille",
    "Heck",
    "Sonogashira",
    "Buchwald-Hartwig",
    "Ullmann-Goldberg",
    "N-arylation",
    "decarboxylative_coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the formation of a bi-heterocyclic system via a specified list of cross-coupling reactions. The function verifies that at least two different heterocycles from a predefined list are present in the reactants and are coupled to form the final product.
    """
    found_strategy = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth=0):
        nonlocal found_strategy, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                if not product_smiles or len(reactants_smiles) < 2:
                    for child in node.get("children", []):
                        dfs_traverse(child, depth)
                    return

                is_coupling = False
                for reaction_type in COUPLING_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_coupling = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

                if is_coupling:
                    heterocycle_in_reactants = {}
                    for i, reactant_smiles in enumerate(reactants_smiles):
                        if not reactant_smiles:
                            continue

                        for heterocycle in HETEROCYCLES_OF_INTEREST:
                            if checker.check_ring(heterocycle, reactant_smiles):
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                                heterocycle_in_reactants[heterocycle] = (
                                    heterocycle_in_reactants.get(heterocycle, []) + [i]
                                )

                    if len(heterocycle_in_reactants) >= 2:
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "distinct_heterocycles_in_reactants_of_coupling_step",
                                "operator": ">=",
                                "value": 2
                            }
                        })
                        heterocycles_in_product = []
                        for heterocycle in heterocycle_in_reactants:
                            if checker.check_ring(heterocycle, product_smiles):
                                # Only add to findings_json if not already present to avoid duplicates
                                if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                                heterocycles_in_product.append(heterocycle)

                        if len(heterocycles_in_product) >= 2:
                            findings_json["structural_constraints"].append({
                                "type": "count",
                                "details": {
                                    "target": "distinct_heterocycles_in_product_of_coupling_step",
                                    "operator": ">=",
                                    "value": 2
                                }
                            })
                            found_strategy = True
            except Exception:
                pass

        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    return found_strategy, findings_json

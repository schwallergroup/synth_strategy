from typing import Tuple, Dict, List
import copy
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

HETEROARYL_COUPLING_REACTIONS = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki",
    "Stille",
    "Negishi",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Ullmann condensation",
    "Buchwald-Hartwig",
    "N-arylation",
    "Sonogashira",
    "decarboxylative_coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage (depth <= 2) coupling reaction between a pyrimidine-containing reactant and a pyridine-containing reactant.
    The product must contain both ring systems. The reaction is identified by checking against a specific list of named coupling reactions.
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

    found_specific_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_specific_coupling, findings_json

        if node.get("type") == "reaction" and depth <= 2:
            # Record positional constraint
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "reaction", "position": "depth <= 2"}})

            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                has_pyrimidine_product = checker.check_ring("pyrimidine", product_smiles)
                if has_pyrimidine_product:
                    if "pyrimidine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("pyrimidine")

                has_pyridine_product = checker.check_ring("pyridine", product_smiles)
                if has_pyridine_product:
                    if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("pyridine")

                if has_pyrimidine_product and has_pyridine_product:
                    # Record co-occurrence in product
                    findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["pyrimidine", "pyridine"], "scope": "product"}})

                    # Check that both pyrimidine and pyridine fragments are present in the reactants
                    pyrimidine_in_reactants = False
                    for r in reactants_smiles:
                        if checker.check_ring("pyrimidine", r):
                            pyrimidine_in_reactants = True
                            if "pyrimidine" not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append("pyrimidine")
                            break

                    pyridine_in_reactants = False
                    for r in reactants_smiles:
                        if checker.check_ring("pyridine", r):
                            pyridine_in_reactants = True
                            if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append("pyridine")
                            break

                    if pyrimidine_in_reactants and pyridine_in_reactants:
                        # Record co-occurrence in reactants
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["pyrimidine", "pyridine"], "scope": "reactants"}})

                        # Check if the reaction is a known coupling type from the predefined list
                        for rxn_type in HETEROARYL_COUPLING_REACTIONS:
                            if checker.check_reaction(rxn_type, rsmi):
                                found_specific_coupling = True
                                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                                return  # Found a match, no need to check further

            except (KeyError, IndexError):
                # Silently ignore errors from malformed reaction data
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node.get("type") != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    return found_specific_coupling, findings_json
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


HETEROCYCLE_RINGS = [
    "furan", "pyran", "pyrrole", "pyridine", "pyrazole", "imidazole",
    "oxazole", "thiazole", "pyrimidine", "pyrazine", "pyridazine", "triazole",
    "tetrazole", "indole", "quinoline", "isoquinoline", "purine", "thiophene",
    "isoxazole", "isothiazole", "oxadiazole", "thiadiazole", "benzoxazole",
    "benzothiazole", "benzimidazole", "indazole", "benzotriazole",
]

COUPLING_REACTIONS = [
    "Suzuki", "Negishi", "Stille", "Buchwald-Hartwig", "N-arylation", "Heck",
    "Sonogashira", "Ullmann-Goldberg", "decarboxylative_coupling",
]

HETEROCYCLE_FORMATION_REACTIONS = [
    "benzimidazole_derivatives_aldehyde",
    "benzimidazole_derivatives_carboxylic-acid/ester", "benzothiazole",
    "benzoxazole_arom-aldehyde", "benzoxazole_carboxylic-acid", "imidazole",
    "thiazole", "tetrazole_terminal", "pyrazole", "oxadiazole",
    "Pictet-Spengler", "Fischer indole", "benzofuran", "benzothiophene",
    "indole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route involves a convergent reaction step (i.e., with two or more reactants) that is either a named heterocycle formation reaction or a named coupling reaction involving a heterocyclic fragment. The specific heterocycles and reaction names are defined in enumerated module-level lists.
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

    heterocycle_coupling_found = False

    # In a full implementation, max_depth would be calculated here.
    max_depth = 10 # Placeholder

    def dfs_traverse(node, depth, max_depth):
        nonlocal heterocycle_coupling_found, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_part, _, product_part = rsmi.partition('>>')
            reactants = reactants_part.split('.')

            # Strategy requires a convergent step (>= 2 reactants)
            if len(reactants) >= 2:
                is_coupling = False
                for name in COUPLING_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        is_coupling = True
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)
                        break

                is_heterocycle_formation = False
                for name in HETEROCYCLE_FORMATION_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        is_heterocycle_formation = True
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)
                        break

                if is_heterocycle_formation:
                    # A convergent named heterocycle formation reaction is a positive hit.
                    heterocycle_coupling_found = True
                    # Record structural constraint for heterocycle formation
                    if {"type": "co-occurrence", "details": {"targets": ["convergent_reaction_step", "named_heterocycle_formation_reaction"], "description": "A single reaction step must be convergent (>= 2 reactants) and be a named heterocycle formation reaction."}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["convergent_reaction_step", "named_heterocycle_formation_reaction"], "description": "A single reaction step must be convergent (>= 2 reactants) and be a named heterocycle formation reaction."}})

                elif is_coupling:
                    # For coupling, ensure at least one reactant is a heterocycle and a heterocycle is in the product.
                    reactant_has_heterocycle = False
                    for r in reactants:
                        for ring in HETEROCYCLE_RINGS:
                            if checker.check_ring(ring, r):
                                reactant_has_heterocycle = True
                                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                                break
                        if reactant_has_heterocycle: # Optimization: stop checking reactants if one heterocycle found
                            break

                    if reactant_has_heterocycle:
                        product_has_heterocycle = False
                        for ring in HETEROCYCLE_RINGS:
                            if checker.check_ring(ring, product_part):
                                product_has_heterocycle = True
                                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                                break

                        if product_has_heterocycle:
                            heterocycle_coupling_found = True
                            # Record structural constraint for coupling with heterocycles
                            if {"type": "co-occurrence", "details": {"targets": ["convergent_reaction_step", "named_coupling_reaction", "heterocycle_in_reactant", "heterocycle_in_product"], "description": "A single reaction step must be a convergent (>= 2 reactants) coupling reaction where at least one reactant contains a heterocycle and the product also contains a heterocycle."}} not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["convergent_reaction_step", "named_coupling_reaction", "heterocycle_in_reactant", "heterocycle_in_product"], "description": "A single reaction step must be a convergent (>= 2 reactants) coupling reaction where at least one reactant contains a heterocycle and the product also contains a heterocycle."}})

        for child in node.get("children", []):
            if heterocycle_coupling_found:
                return # Optimization: stop traversal once found
            
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth, max_depth)

    dfs_traverse(route, 0, max_depth)
    return heterocycle_coupling_found, findings_json

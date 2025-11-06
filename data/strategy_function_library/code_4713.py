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


# Refactoring for Enumeration: Isolate lists as module-level constants
NITROGEN_HETEROCYCLES = [
    "pyridine", "pyrrole", "imidazole", "pyrazole", "triazole", "tetrazole",
    "pyrimidine", "pyrazine", "indole", "quinoline", "isoquinoline", "purine",
    "benzimidazole", "benzotriazole", "piperidine", "piperazine", "morpholine",
]

COUPLING_REACTIONS = [
    "Suzuki", "Buchwald-Hartwig", "Stille", "Negishi", "Heck", "Sonogashira",
    "N-arylation", "Ullmann-Goldberg", "Ullmann condensation", "Chan-Lam",
    "Goldberg coupling", "decarboxylative_coupling",
]

HETEROCYCLE_FORMATION_REACTIONS = [
    "benzimidazole_derivatives_carboxylic-acid/ester", "benzimidazole_derivatives_aldehyde",
    "benzothiazole", "benzoxazole_arom-aldehyde", "benzoxazole_carboxylic-acid",
    "thiazole", "tetrazole_terminal", "tetrazole_connect_regioisomere_1",
    "tetrazole_connect_regioisomere_2", "1,2,4-triazole_acetohydrazide",
    "1,2,4-triazole_carboxylic-acid/ester", "pyrazole", "Paal-Knorr pyrrole",
    "triaryl-imidazole", "Fischer indole", "imidazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects synthesis steps that either form a new nitrogen-containing heterocycle
    via a named reaction or couple fragments using a named coupling reaction, where
    at least one fragment contains a nitrogen heterocycle. The specific heterocycles
    and reaction names are defined in the module-level lists: NITROGEN_HETEROCYCLES,
    COUPLING_REACTIONS, and HETEROCYCLE_FORMATION_REACTIONS.
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

    # Track coupling reactions that join heterocycles
    heterocycle_couplings = []

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_couplings, findings_json
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                is_coupling = False
                coupling_type = None
                for rxn_name in COUPLING_REACTIONS:
                    if checker.check_reaction(rxn_name, rsmi):
                        is_coupling = True
                        coupling_type = rxn_name
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        break

                is_heterocycle_formation = False
                for rxn_name in HETEROCYCLE_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn_name, rsmi):
                        is_heterocycle_formation = True
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        break

                reactant_heterocycles = []
                for reactant in reactants:
                    reactant_rings = []
                    for ring_name in NITROGEN_HETEROCYCLES:
                        if checker.check_ring(ring_name, reactant):
                            reactant_rings.append(ring_name)
                            if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                    if reactant_rings:
                        reactant_heterocycles.append((reactant, reactant_rings))

                product_rings = []
                for ring_name in NITROGEN_HETEROCYCLES:
                    if checker.check_ring(ring_name, product):
                        product_rings.append(ring_name)
                        if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring_name)

                # Case 1: Coupling reaction joining fragments with heterocycles
                if is_coupling and len(reactant_heterocycles) >= 1 and product_rings:
                    max_reactant_rings = (
                        max([len(rings) for _, rings in reactant_heterocycles])
                        if reactant_heterocycles
                        else 0
                    )
                    if len(reactant_heterocycles) >= 2 or len(product_rings) > max_reactant_rings:
                        heterocycle_couplings.append(
                            (reactant_heterocycles, product_rings, coupling_type)
                        )
                        # Add structural constraint for coupling
                        if {"type": "co-occurrence", "details": {"targets": ["coupling_reaction", "nitrogen_heterocycle"], "description": "A named coupling reaction occurs where at least one reactant contains a nitrogen heterocycle."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["coupling_reaction", "nitrogen_heterocycle"], "description": "A named coupling reaction occurs where at least one reactant contains a nitrogen heterocycle."}})

                # Case 2: Heterocycle formation reaction that creates a new heterocycle
                elif is_heterocycle_formation and product_rings:
                    reactant_ring_set = set()
                    for _, rings in reactant_heterocycles:
                        reactant_ring_set.update(rings)
                    new_rings = [ring for ring in product_rings if ring not in reactant_ring_set]
                    if new_rings:
                        heterocycle_couplings.append(
                            (reactant_heterocycles, product_rings, "formation")
                        )
                        # Add structural constraint for heterocycle formation
                        if {"type": "co-occurrence", "details": {"targets": ["heterocycle_formation_reaction", "ring_formation"], "description": "A named heterocycle formation reaction results in the creation of a new nitrogen heterocycle not present in the reactants."}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["heterocycle_formation_reaction", "ring_formation"], "description": "A named heterocycle formation reaction results in the creation of a new nitrogen heterocycle not present in the reactants."}})
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

            except Exception:
                # Silently ignore errors in reaction analysis
                pass

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other non-reaction type
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = len(heterocycle_couplings) > 0

    # Add the final count constraint if the result is True
    if result:
        if {"type": "count", "details": {"target": "heterocycle_coupling_or_formation_event", "operator": ">", "value": 0, "description": "The overall strategy is valid if at least one of the defined co-occurrence events is found in the synthesis route."}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "heterocycle_coupling_or_formation_event", "operator": ">", "value": 0, "description": "The overall strategy is valid if at least one of the defined co-occurrence events is found in the synthesis route."}})

    return result, findings_json

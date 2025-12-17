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


# Refactoring for Enumeration: Isolate the lists of reactions and rings.
HETEROCYCLE_FORMATION_REACTIONS = [
    "Formation of NOS Heterocycles",
    "Paal-Knorr pyrrole synthesis",
    "benzimidazole_derivatives_carboxylic-acid/ester",
    "benzimidazole_derivatives_aldehyde",
    "benzothiazole",
    "benzoxazole_arom-aldehyde",
    "benzoxazole_carboxylic-acid",
    "thiazole",
    "tetrazole_terminal",
    "tetrazole_connect_regioisomere_1",
    "tetrazole_connect_regioisomere_2",
    "1,2,4-triazole_acetohydrazide",
    "1,2,4-triazole_carboxylic-acid/ester",
    "pyrazole",
    "Fischer indole",
    "benzofuran",
    "benzothiophene",
    "indole",
    "oxadiazole",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Pyrazole formation",
]

HETEROCYCLES_OF_INTEREST = [
    "furan", "pyrrole", "thiophene", "pyridine", "oxazole", "thiazole",
    "imidazole", "pyrimidine", "pyrazine", "pyridazine", "triazole",
    "tetrazole", "oxadiazole", "thiadiazole", "isoxazole", "isothiazole",
    "benzoxazole", "benzothiazole", "benzimidazole", "indole", "quinoline",
    "isoquinoline", "purine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route contains late-stage heterocycle formation.
    Late stage is defined as occurring in the first half of the synthesis (low depth).
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

    max_depth = 0
    heterocycle_formation_depths = []
    result = False

    # First pass to find max depth
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            # Determine the depth for the recursive call based on node type
            if node["type"] == "reaction":
                # If current node is 'reaction', depth remains the same for children
                next_depth = current_depth
            else:
                # If current node is 'chemical', depth increases for children
                next_depth = current_depth + 1
            find_max_depth(child, next_depth)

    # Second pass to find heterocycle formations
    def find_heterocycle_formation(node, depth=0):
        nonlocal result
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            is_heterocycle_formation_reaction = False
            for rxn_type in HETEROCYCLE_FORMATION_REACTIONS:
                if checker.check_reaction(rxn_type, rsmi):
                    is_heterocycle_formation_reaction = True
                    if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    break

            if not is_heterocycle_formation_reaction:
                try:
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    product_mol = Chem.MolFromSmiles(product) if product else None

                    if product_mol and all(reactant_mols):
                        for ring_name in HETEROCYCLES_OF_INTEREST:
                            if checker.check_ring(ring_name, product):
                                if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                                reactant_has_ring = False
                                for r in reactants:
                                    if checker.check_ring(ring_name, r):
                                        reactant_has_ring = True
                                        break

                                if not reactant_has_ring:
                                    is_heterocycle_formation_reaction = True
                                    # Add a generic 'ring_formation' reaction if not already present
                                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                                    break
                except Exception:
                    pass

            if is_heterocycle_formation_reaction:
                heterocycle_formation_depths.append(depth)

        for child in node.get("children", []):
            # Determine the depth for the recursive call based on node type
            if node["type"] == "reaction":
                # If current node is 'reaction', depth remains the same for children
                next_depth = depth
            else:
                # If current node is 'chemical', depth increases for children
                next_depth = depth + 1
            find_heterocycle_formation(child, next_depth)

    find_max_depth(route)
    find_heterocycle_formation(route)

    if heterocycle_formation_depths:
        heterocycle_formation_depths.sort()
        earliest_formation_depth = heterocycle_formation_depths[0]

        is_late_stage = earliest_formation_depth <= max_depth / 2
        result = is_late_stage

        if result:
            # Add the structural constraint if the condition is met
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "heterocycle_formation",
                    "position": "in_latter_half_of_route",
                    "condition": "The earliest heterocycle formation must occur at a depth less than or equal to half the maximum route depth (depth <= max_depth / 2)."
                }
            })

    return result, findings_json

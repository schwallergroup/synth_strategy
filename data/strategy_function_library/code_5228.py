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


HETEROCYCLE_TYPES = [
    "pyridine", "pyrimidine", "pyrazine", "pyridazine", "triazine", "tetrazine",
    "pyrrole", "furan", "thiophene", "imidazole", "oxazole", "thiazole",
    "pyrazole", "isoxazole", "isothiazole", "triazole", "tetrazole", "indole",
    "benzofuran", "benzothiophene", "benzimidazole", "benzoxazole",
    "benzothiazole", "purine", "quinoline", "isoquinoline", "phthalazine",
    "quinoxaline", "quinazoline", "cinnoline",
]

FUNCTIONALIZATION_REACTIONS = [
    "Acylation", "Alkylation", "N-alkylation", "O-alkylation", "S-alkylation",
    "Amination", "Buchwald-Hartwig", "Suzuki", "Sonogashira", "Heck", "Negishi",
    "Stille", "Methylation", "Friedel-Crafts", "Halogenation", "Nitration",
    "Sulfonation", "N-arylation", "C-alkylation", "C-arylation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis follows a linear build-up approach where a core heterocycle
    is sequentially functionalized.
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

    # Track reactions and their depths
    reactions_by_depth = {}

    # Store information about the core heterocycle
    core_info = {
        "found": False,
        "type": None,
        "molecules": [],  # Track all molecules containing the heterocycle
        "target_mol": None,  # The final target molecule
    }

    # Track functionalization steps
    functionalization_steps = []

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # If we haven't identified a heterocycle yet, check this molecule
            if not core_info["found"]:
                for heterocycle in HETEROCYCLE_TYPES:
                    if checker.check_ring(heterocycle, mol_smiles):
                        core_info["found"] = True
                        core_info["type"] = heterocycle
                        core_info["molecules"].append((depth, mol_smiles))
                        if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                        # The first molecule we encounter is the target
                        if core_info["target_mol"] is None:
                            core_info["target_mol"] = mol_smiles
                        break

            # If we already found a heterocycle, check if this molecule also contains it
            elif core_info["found"]:
                if checker.check_ring(core_info["type"], mol_smiles):
                    core_info["molecules"].append((depth, mol_smiles))
                    if core_info["type"] not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(core_info["type"])

        elif node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactions_by_depth[depth] = rsmi

                # Extract product and reactants
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this reaction involves the heterocycle
                if core_info["found"]:
                    product_has_heterocycle = checker.check_ring(core_info["type"], product)

                    # Check if at least one reactant has the heterocycle
                    reactant_with_heterocycle = False
                    for reactant in reactants:
                        if checker.check_ring(core_info["type"], reactant):
                            reactant_with_heterocycle = True
                            break

                    # If both product and at least one reactant have the heterocycle,
                    # this might be a functionalization reaction
                    if product_has_heterocycle and reactant_with_heterocycle:
                        # Record the co-occurrence constraint
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "core_heterocycle_in_product",
                                    "core_heterocycle_in_reactant"
                                ],
                                "scope": "per_reaction",
                                "description": "A reaction is only considered a valid functionalization step if the identified core heterocycle is present in both a reactant and the product."
                            }
                        })
                        # Check if this is a known functionalization reaction type
                        for rxn_type in FUNCTIONALIZATION_REACTIONS:
                            if checker.check_reaction(rxn_type, rsmi):
                                functionalization_steps.append((depth, rxn_type, rsmi))
                                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                                break

        for child in node.get("children", []):
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start the traversal
    dfs_traverse(route)

    # Sort molecules by depth (ascending)
    core_info["molecules"].sort(key=lambda x: x[0])

    # Check if we have a linear sequence with at least 2 functionalization steps
    is_linear = len(functionalization_steps) >= 2
    if is_linear:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "functionalization_reaction",
                "operator": ">=",
                "value": 2,
                "description": "The route must contain at least two valid functionalization reactions."
            }
        })

    # Check if we found and maintained a heterocycle through multiple steps
    heterocycle_maintained = core_info["found"] and len(core_info["molecules"]) >= 3
    if heterocycle_maintained:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "heterocycle_bearing_molecule",
                "operator": ">=",
                "value": 3,
                "description": "The identified core heterocycle must be present in at least three molecules along the main synthetic path."
            }
        })

    # Check if the reactions form a linear sequence
    sequential_modification = True
    if is_linear and heterocycle_maintained:
        # Sort functionalization steps by depth
        functionalization_steps.sort(key=lambda x: x[0])

        # In retrosynthetic analysis, we need to check if the sequence is continuous
        # Each step should modify the heterocycle from the previous step
        for i in range(1, len(functionalization_steps)):
            current_depth = functionalization_steps[i][0]
            prev_depth = functionalization_steps[i - 1][0]

            # Check if there are no gaps in the sequence
            # In a linear sequence, the depths should be consecutive or have a small gap
            if current_depth - prev_depth > 2:  # Allow for some flexibility
                sequential_modification = False
                break
        if sequential_modification:
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "events": [
                        "functionalization_reaction",
                        "functionalization_reaction"
                    ],
                    "max_depth_difference": 2,
                    "description": "Consecutive functionalization steps must be separated by at most 2 levels in the synthesis tree."
                }
            })

    result = is_linear and heterocycle_maintained and sequential_modification

    return result, findings_json

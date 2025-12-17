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


COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki",
    "Negishi",
    "Stille",
    "Heck",
    "Sonogashira",
    "Buchwald-Hartwig",
    "N-arylation",
    "Ullmann-Goldberg",
    "Ullmann condensation",
    "Chan-Lam",
    "Kumada",
    "Hiyama-Denmark",
    "decarboxylative_coupling",
    "Michael addition",
    "Minisci",
    "Friedel-Crafts alkylation",
    "Catellani reaction",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy involving late-stage coupling with a heterocyclic fragment
    (specifically barbituric acid or similar structures).
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

    found_heterocycle_coupling = False

    # List of common heterocycles to check
    heterocycles = [
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazine",
        "furan",
        "thiophene",
        "pyrrole",
        "imidazole",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "indole",
        "benzimidazole",
        "purine",
        "quinoline",
        "isoquinoline",
        "benzoxazole",
        "benzothiazole",
        "morpholine",
        "piperidine",
        "piperazine",
        "pyrrolidine",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycle_coupling, findings_json

        if node["type"] == "reaction" and depth <= 2:  # Final, penultimate, or antepenultimate step
            # Structural constraint: positional - reaction_of_interest, late_stage
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "reaction_of_interest",
                    "position": "late_stage",
                    "description": "The key reaction must occur within the last three steps (depth <= 2)."
                }
            })
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                if len(reactants_smiles) > 1:
                    heterocycle_reactants = []
                    heterocycle_types = []

                    # Check specifically for barbituric acid or derivatives
                    barbituric_found = False
                    for idx, reactant_smi in enumerate(reactants_smiles):
                        # Check for pyrimidine ring with multiple amides/urea groups
                        if checker.check_ring("pyrimidine", reactant_smi):
                            if "pyrimidine" not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append("pyrimidine")

                            urea_found = checker.check_fg("Urea", reactant_smi)
                            amide_count = sum(1 for _ in range(3) if checker.check_fg("Amide", reactant_smi))

                            if urea_found:
                                if "Urea" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Urea")
                            if amide_count >= 1:
                                if "Amide" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Amide")

                            if urea_found or (amide_count >= 2):
                                barbituric_found = True
                                if idx not in heterocycle_reactants:
                                    heterocycle_reactants.append(idx)
                                    heterocycle_types.append("pyrimidine")
                                # Structural constraint: co-occurrence - reactant (pyrimidine + Urea/Amide)
                                findings_json["structural_constraints"].append({
                                    "type": "co-occurrence",
                                    "details": {
                                        "scope": "reactant",
                                        "description": "A key reactant must contain a pyrimidine ring and also either a Urea group or at least two Amide groups."
                                    }
                                })
                                if amide_count >= 2:
                                    # Structural constraint: count - Amide
                                    findings_json["structural_constraints"].append({
                                        "type": "count",
                                        "details": {
                                            "target": "Amide",
                                            "scope": "reactant",
                                            "operator": ">=",
                                            "value": 2
                                        }
                                    })

                    if len(heterocycle_reactants) >= 2 or (
                        len(heterocycle_reactants) >= 1 and barbituric_found
                    ):
                        is_coupling = False
                        for rxn_type in COUPLING_REACTIONS_OF_INTEREST:
                            if checker.check_reaction(rxn_type, rsmi):
                                is_coupling = True
                                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                                break

                        if is_coupling:
                            heterocycles_in_product = sum(
                                1
                                for ht in heterocycle_types
                                if checker.check_ring(ht, product_smiles)
                            )
                            if heterocycles_in_product >= len(heterocycle_types):
                                found_heterocycle_coupling = True
                                # Structural constraint: structural - ring_conservation
                                findings_json["structural_constraints"].append({
                                    "type": "structural",
                                    "details": {
                                        "scope": "reaction",
                                        "target": "pyrimidine",
                                        "constraint": "ring_conservation",
                                        "description": "The pyrimidine ring from the key reactant must be conserved in the product."
                                    }
                                })
                                # Structural constraint: co-occurrence - reaction (coupling)
                                findings_json["structural_constraints"].append({
                                    "type": "co-occurrence",
                                    "details": {
                                        "scope": "reaction",
                                        "description": "The key reaction must be either one of the specified coupling reactions, or a non-coupling reaction that forms a single product molecule."
                                    }
                                })
                        else:
                            heterocycles_in_product = sum(
                                1
                                for ht in heterocycle_types
                                if checker.check_ring(ht, product_smiles)
                            )

                            if heterocycles_in_product >= len(heterocycle_types):
                                if "." not in product_smiles:
                                    found_heterocycle_coupling = True
                                    # Structural constraint: structural - ring_conservation
                                    findings_json["structural_constraints"].append({
                                        "type": "structural",
                                        "details": {
                                            "scope": "reaction",
                                            "target": "pyrimidine",
                                            "constraint": "ring_conservation",
                                            "description": "The pyrimidine ring from the key reactant must be conserved in the product."
                                        }
                                    })
                                    # Structural constraint: co-occurrence - reaction (non-coupling single product)
                                    findings_json["structural_constraints"].append({
                                        "type": "co-occurrence",
                                        "details": {
                                            "scope": "reaction",
                                            "description": "The key reaction must be either one of the specified coupling reactions, or a non-coupling reaction that forms a single product molecule."
                                        }
                                    })
            except Exception:
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    return found_heterocycle_coupling, findings_json

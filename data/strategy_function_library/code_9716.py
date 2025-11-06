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


HETEROCYCLE_FORMATION_DATA = [
    (
        "benzimidazole",
        [
            "benzimidazole_derivatives_aldehyde",
            "benzimidazole_derivatives_carboxylic-acid/ester",
        ],
    ),
    ("benzoxazole", ["benzoxazole_arom-aldehyde", "benzoxazole_carboxylic-acid"]),
    ("benzothiazole", ["benzothiazole"]),
    ("indole", ["Fischer indole", "indole"]),
    ("imidazole", ["imidazole", "triaryl-imidazole"]),
    ("pyrrole", ["Paal-Knorr pyrrole", "pyrrole"]),
    (
        "tetrazole",
        [
            "tetrazole_terminal",
            "tetrazole_connect_regioisomere_1",
            "tetrazole_connect_regioisomere_2",
        ],
    ),
    ("thiazole", ["thiazole"]),
    ("pyrazole", ["pyrazole"]),
    ("oxazole", ["oxazole"]),
    (
        "triazole",
        ["triazole", "Huisgen_Cu-catalyzed_1,4-subst", "Huisgen_Ru-catalyzed_1,5_subst"],
    ),
    ("furan", ["furan", "benzofuran"]),
    ("thiophene", ["thiophene", "benzothiophene"]),
    ("pyridine", ["pyridine", "3-nitrile-pyridine"]),
    ("quinoline", ["quinoline", "Friedlaender chinoline"]),
    ("isoquinoline", ["isoquinoline"]),
    ("oxadiazole", ["oxadiazole", "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole"]),
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage (final or penultimate step) formation of a specific heterocycle from precursors that do not contain that heterocycle.
    The formation must be identified as one of the specified named reactions associated with the target heterocycle, as defined in HETEROCYCLE_FORMATION_DATA.
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

    # Track if we found the heterocycle formation
    heterocycle_formed = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formed, findings_json

        if node["type"] == "reaction" and depth < 2:
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for heterocycle formation via a specific named reaction
            for heterocycle, reaction_types in HETEROCYCLE_FORMATION_DATA:
                if checker.check_ring(heterocycle, product_smiles):
                    if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                    reactants_have_heterocycle = any(
                        checker.check_ring(heterocycle, r) for r in reactants_smiles
                    )

                    if not reactants_have_heterocycle:
                        # This is a confirmed formation. Now, check if it's via a known reaction type.
                        # Add negation constraint
                        if {"type": "negation", "details": {"target": "target_heterocycle_in_reactants", "scope": "reaction_step"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "target_heterocycle_in_reactants", "scope": "reaction_step"}})

                        for rxn_type in reaction_types:
                            if checker.check_reaction(rxn_type, rsmi):
                                heterocycle_formed = True
                                if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                                # Add co-occurrence constraint
                                if {"type": "co-occurrence", "details": {"targets": ["target_heterocycle_in_product", "named_heterocycle_formation_reaction"], "scope": "reaction_step"}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["target_heterocycle_in_product", "named_heterocycle_formation_reaction"], "scope": "reaction_step"}})
                                # Add positional constraint if applicable
                                if depth == 0 or depth == 1: # Final or penultimate step
                                    if {"type": "positional", "details": {"target": "heterocycle_formation_reaction", "position": "last_or_penultimate_stage"}} not in findings_json["structural_constraints"]:
                                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "heterocycle_formation_reaction", "position": "last_or_penultimate_stage"}})

                                # Found a valid reaction, no need to check other types for this heterocycle
                                break
                
                if heterocycle_formed:
                    # Found a match, no need to check other heterocycles in this reaction
                    break

        # Traverse children
        for child in node.get("children", []):
            if heterocycle_formed: # Optimization: stop traversing if already found
                return
            
            # New depth calculation logic
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other non-reaction type
                dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return heterocycle_formed, findings_json

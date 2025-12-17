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
    "pyrrole", "pyridine", "pyrazole", "imidazole", "oxazole", "thiazole",
    "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole", "indole",
    "quinoline", "isoquinoline", "benzoxazole", "benzothiazole", "benzimidazole",
    "furan", "thiophene", "oxadiazole", "thiadiazole", "isoxazole",
    "quinazoline", "pteridin", "purine",
]

HETEROCYCLE_FORMATION_REACTIONS = [
    "Paal-Knorr pyrrole synthesis", "Fischer indole", "benzimidazole_derivatives_aldehyde",
    "benzothiazole", "benzoxazole_arom-aldehyde", "benzoxazole_carboxylic-acid",
    "thiazole", "tetrazole_terminal", "Huisgen_Cu-catalyzed_1,4-subst",
    "Huisgen_Ru-catalyzed_1,5_subst", "1,2,4-triazole_acetohydrazide", "pyrazole",
    "oxadiazole", "Pictet-Spengler", "benzimidazole_derivatives_carboxylic-acid/ester",
    "Formation of NOS Heterocycles", "Niementowski_quinazoline", "imidazole",
    "triaryl-imidazole", "benzofuran", "benzothiophene", "indole",
    "tetrazole_connect_regioisomere_1", "tetrazole_connect_regioisomere_2",
    "Huisgen_disubst-alkyne", "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition", "Huisgen alkene-azide 1,3 dipolar cycloaddition",
]

FUNCTIONAL_GROUPS_OF_INTEREST = [
    "Nitro group", "Nitrile", "Primary amine", "Secondary amine", "Tertiary amine",
    "Carboxylic acid", "Ester", "Amide", "Aldehyde", "Ketone", "Alcohol",
    "Aromatic halide", "Primary halide", "Secondary halide", "Tertiary halide",
    "Azide", "Acyl halide", "Anhydride", "Sulfonamide", "Sulfone", "Sulfoxide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy involving late-stage formation of specific heterocycles, as defined in `HETEROCYCLE_FORMATION_REACTIONS` or `HETEROCYCLE_RINGS`. This strategy requires that the synthesis preserves an aromatic core and includes modifications of functional groups listed in `FUNCTIONAL_GROUPS_OF_INTEREST` during the early stages.
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

    final_step_heterocycle_formation = False
    early_fg_modifications = False
    aromatic_core_preserved = True

    def dfs_traverse(node, depth=0):
        nonlocal final_step_heterocycle_formation, early_fg_modifications, aromatic_core_preserved, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]
            reactants = reactants_part.split(".")

            heterocycle_formed = False
            for reaction_type in HETEROCYCLE_FORMATION_REACTIONS:
                if checker.check_reaction(reaction_type, rsmi):
                    heterocycle_formed = True
                    if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                    if "heterocycle_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("heterocycle_formation")
                    break

            if not heterocycle_formed:
                product_has_heterocycle = False
                product_heterocycles = []
                for ring in HETEROCYCLE_RINGS:
                    if checker.check_ring(ring, product_part):
                        product_has_heterocycle = True
                        product_heterocycles.append(ring)
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)

                if product_has_heterocycle:
                    for ring in product_heterocycles:
                        ring_in_reactants = False
                        for reactant in reactants:
                            if checker.check_ring(ring, reactant):
                                ring_in_reactants = True
                                break
                        if not ring_in_reactants:
                            heterocycle_formed = True
                            if "heterocycle_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("heterocycle_formation")
                            break

            fg_modified = False
            for fg in FUNCTIONAL_GROUPS_OF_INTEREST:
                reactant_has_fg = any(checker.check_fg(fg, r) for r in reactants)
                product_has_fg = checker.check_fg(fg, product_part)
                if reactant_has_fg != product_has_fg:
                    fg_modified = True
                    if fg not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(fg)
                    if "functional_group_modification" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("functional_group_modification")
                    break

            product_has_aromatic = checker.check_aromaticity(product_part)
            reactants_have_aromatic = any(checker.check_aromaticity(r) for r in reactants)
            if reactants_have_aromatic and not product_has_aromatic:
                aromatic_core_preserved = False
                if "aromatic_ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("aromatic_ring_destruction")
                if "aromatic_ring" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("aromatic_ring")

            if depth <= 1:
                if heterocycle_formed:
                    final_step_heterocycle_formation = True
                    # Add structural constraint for late-stage heterocycle formation
                    constraint = {"type": "positional", "details": {"target": "heterocycle_formation", "position": "late_stage"}}
                    if constraint not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(constraint)
            else:
                if fg_modified:
                    early_fg_modifications = True
                    # Add structural constraint for early-stage functional group modification
                    constraint = {"type": "positional", "details": {"target": "functional_group_modification", "position": "early_stage"}}
                    if constraint not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(constraint)

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = final_step_heterocycle_formation and early_fg_modifications and aromatic_core_preserved

    # Add structural constraint for aromatic core preservation (negation)
    if not aromatic_core_preserved:
        constraint = {"type": "negation", "details": {"target": "aromatic_ring_destruction"}}
        # This constraint is met if aromatic_core_preserved is FALSE, meaning destruction occurred.
        # The strategy requires aromatic_core_preserved to be TRUE, so if it's FALSE, this constraint is violated.
        # We only add it to findings_json if the condition for the overall strategy is met (i.e., aromatic_core_preserved is TRUE)
        # However, the instruction is to add 'only the elements that were actually detected'.
        # If aromatic_core_preserved is FALSE, then 'aromatic_ring_destruction' was detected.
        # The structural constraint is about the *absence* of destruction. So if destruction happens, this constraint is *not* met.
        # The strategy is looking for the *preservation* of the aromatic core, which means the *absence* of aromatic_ring_destruction.
        # So, if aromatic_core_preserved is TRUE, then the negation constraint is effectively met.
        # Let's re-evaluate: The strategy says 'requires that the synthesis preserves an aromatic core'.
        # This means 'aromatic_ring_destruction' should NOT happen. So if 'aromatic_ring_destruction' is detected, the strategy fails.
        # The structural constraint in the JSON is 'negation' of 'aromatic_ring_destruction'.
        # If 'aromatic_core_preserved' is TRUE, it means 'aromatic_ring_destruction' did NOT happen.
        # So, if 'aromatic_core_preserved' is TRUE, we should add the negation constraint to findings_json.
        pass # This is handled by the overall result. The atomic check 'aromatic_ring_destruction' is added if it occurs.

    # The structural constraint for 'negation' of 'aromatic_ring_destruction' is met if aromatic_core_preserved is TRUE.
    # So, if aromatic_core_preserved is TRUE, we add the constraint to findings_json.
    if aromatic_core_preserved:
        constraint = {"type": "negation", "details": {"target": "aromatic_ring_destruction"}}
        if constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(constraint)

    return result, findings_json

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


# Refactored lists for identifying late-stage diversification
COMPLEX_FGS_LSD = [
    "Triflate",
    "Tosylate",
    "Mesylate",
    "Sulfonamide",
    "Trifluoro group",
    "Trichloro group",
    "Carbamic ester",
    "Carbamic acid",
    "Carbonic Ester",
    "Substituted dicarboximide",
    "Boronic ester",
    "Silyl protective group",
    "TMS ether protective group",
    "Phosphate ester",
    "Boronic acid",
    "Aromatic halide",
    "Primary halide",
    "Secondary halide",
    "Tertiary halide",
]

COMPLEX_RINGS_LSD = [
    "morpholine",
    "thiomorpholine",
    "piperazine",
    "piperidine",
    "pyrrolidine",
    "indole",
    "quinoline",
    "isoquinoline",
    "benzothiazole",
    "benzoxazole",
    "benzimidazole",
    "tetrazole",
    "triazole",
    "oxadiazole",
    "thiadiazole",
    "purine",
    "pteridin",
    "porphyrin",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
]

LATE_STAGE_REACTIONS_LSD = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Buchwald-Hartwig",
    "N-arylation",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl halide",
    "Heck terminal vinyl",
    "Heck_terminal_vinyl",
    "Heck_non-terminal_vinyl",
    "Stille reaction_aryl",
    "Stille reaction_vinyl",
    "Stille",
    "Chan-Lam alcohol",
    "Chan-Lam amine",
    "Chan-Lam etherification",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Huisgen_Cu-catalyzed_1,4-subst",
    "Negishi coupling",
    "Negishi",
    "Kumada cross-coupling",
    "Hiyama-Denmark Coupling",
    "Ullmann-Goldberg Substitution amine",
    "Ullmann-Goldberg Substitution thiol",
    "Ullmann-Goldberg Substitution aryl alcohol",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route uses a late-stage diversification
    strategy, where complex substituents are introduced in the final steps.
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

    late_diversification = False

    def dfs_traverse(node, depth=0):
        nonlocal late_diversification, findings_json

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and depth <= 3:  # Late stage (depth 0-3)
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    print("No reaction SMILES found in metadata")
                    return

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                if not product or not reactants:
                    print("Empty product or reactants, skipping")
                    return

                product_complex_fgs = []
                product_complex_rings = []
                reactants_complex_fgs = set()
                reactants_complex_rings = set()

                for fg in COMPLEX_FGS_LSD:
                    if checker.check_fg(fg, product):
                        product_complex_fgs.append(fg)
                        if fg not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(fg)
                        print(f"Found complex FG in product: {fg}")

                for ring in COMPLEX_RINGS_LSD:
                    if checker.check_ring(ring, product):
                        product_complex_rings.append(ring)
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                        print(f"Found complex ring in product: {ring}")

                for reactant in reactants:
                    for fg in COMPLEX_FGS_LSD:
                        if checker.check_fg(fg, reactant):
                            reactants_complex_fgs.add(fg)
                            if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append(fg)
                            print(f"Found complex FG in reactant: {fg}")

                    for ring in COMPLEX_RINGS_LSD:
                        if checker.check_ring(ring, reactant):
                            reactants_complex_rings.add(ring)
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)
                            print(f"Found complex ring in reactant: {ring}")

                new_fgs = [fg for fg in product_complex_fgs if fg not in reactants_complex_fgs]
                if new_fgs:
                    print(f"New functional groups added: {new_fgs}")
                    # Add structural constraint for new FG formation
                    if {"type": "positional", "details": {"target": "Formation of a new functional group from the COMPLEX_FGS_LSD list", "position": "depth <= 3"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Formation of a new functional group from the COMPLEX_FGS_LSD list", "position": "depth <= 3"}})

                new_rings = [
                    ring for ring in product_complex_rings if ring not in reactants_complex_rings
                ]
                if new_rings:
                    print(f"New ring structures added: {new_rings}")
                    # Add structural constraint for new ring formation
                    if {"type": "positional", "details": {"target": "Formation of a new ring system from the COMPLEX_RINGS_LSD list", "position": "depth <= 3"}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "Formation of a new ring system from the COMPLEX_RINGS_LSD list", "position": "depth <= 3"}})

                reaction_is_late_stage = False
                for rxn in LATE_STAGE_REACTIONS_LSD:
                    if checker.check_reaction(rxn, rsmi):
                        print(f"Late-stage reaction detected: {rxn}")
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        reaction_is_late_stage = True
                        # Add structural constraint for late-stage reaction
                        if {"type": "positional", "details": {"target": "A reaction from the LATE_STAGE_REACTIONS_LSD list", "position": "depth <= 3"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "A reaction from the LATE_STAGE_REACTIONS_LSD list", "position": "depth <= 3"}})
                        break

                if new_fgs or new_rings or reaction_is_late_stage:
                    print(f"Detected late-stage diversification at depth {depth}")
                    late_diversification = True

            except Exception as e:
                print(f"Error in diversification detection: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node['type'] == 'reaction':
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {late_diversification}")
    return late_diversification, findings_json

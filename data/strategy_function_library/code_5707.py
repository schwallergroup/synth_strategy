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


HETEROCYCLES_OF_INTEREST = [
    "furan", "pyran", "pyrrole", "pyridine", "pyrazole", "imidazole",
    "oxazole", "thiazole", "isoxazole", "pyrimidine", "pyrazine", "triazole",
    "tetrazole", "indole", "quinoline", "isoquinoline", "thiophene",
    "benzoxazole", "benzothiazole", "benzimidazole"
]

COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki coupling with boronic acids", "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf", "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters", "Negishi coupling", "Stille reaction_aryl",
    "Stille reaction_vinyl", "Stille reaction_benzyl", "Stille reaction_allyl",
    "Stille reaction_aryl OTf", "Stille reaction_vinyl OTf", "Stille reaction_benzyl OTf",
    "Stille reaction_allyl OTf", "Stille reaction_other", "Stille reaction_other OTf",
    "Heck terminal vinyl", "Heck_terminal_vinyl", "Heck_non-terminal_vinyl",
    "Oxidative Heck reaction", "Oxidative Heck reaction with vinyl ester",
    "Sonogashira acetylene_aryl halide", "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl OTf", "Sonogashira alkyne_aryl OTf",
    "Sonogashira acetylene_alkenyl halide", "Sonogashira alkyne_alkenyl halide",
    "Sonogashira acetylene_acyl halide", "Sonogashira alkyne_acyl halide",
    "Buchwald-Hartwig", "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Kumada cross-coupling", "Hiyama-Denmark Coupling", "Ullmann condensation",
    "Ullmann-Goldberg Substitution amine", "Ullmann-Goldberg Substitution thiol",
    "Ullmann-Goldberg Substitution aryl alcohol", "Goldberg coupling aryl amine-aryl chloride",
    "Goldberg coupling aryl amide-aryl chloride", "Goldberg coupling", "Chan-Lam alcohol",
    "Chan-Lam amine", "Chan-Lam etherification", "Aryllithium cross-coupling",
    "decarboxylative_coupling", "Catellani reaction ortho", "Catellani reaction para"
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy where a specific heterocycle is formed early in the synthesis, followed by a late-stage coupling reaction on the heterocyclic core. The specific heterocycles and coupling reactions are defined in the HETEROCYCLES_OF_INTEREST and COUPLING_REACTIONS_OF_INTEREST lists, respectively.
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

    result = False

    heterocycle_formations = []
    coupling_reactions = []

    final_product_smiles = route["smiles"]

    final_product_heterocycles = []
    for heterocycle in HETEROCYCLES_OF_INTEREST:
        if checker.check_ring(heterocycle, final_product_smiles):
            final_product_heterocycles.append(heterocycle)
            if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

    if not final_product_heterocycles:
        return False, findings_json
    else:
        # Structural constraint: any_heterocycle_of_interest must be present in the final product
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "any_heterocycle_of_interest"
                ],
                "context": "must be present in the final product"
            }
        })

    def traverse(node, depth=0):
        nonlocal heterocycle_formations, coupling_reactions, findings_json
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                product_smiles = rsmi.split(">")[-1]
                reactants_smiles = rsmi.split(">")[0].split(".")

                heterocycle_formed_in_this_step = False
                for heterocycle in HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(heterocycle, product_smiles):
                        reactants_have_heterocycle = False
                        for reactant in reactants_smiles:
                            if checker.check_ring(heterocycle, reactant):
                                reactants_have_heterocycle = True
                                break
                        if not reactants_have_heterocycle:
                            heterocycle_formations.append((depth, heterocycle))
                            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                            if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                            heterocycle_formed_in_this_step = True
                            # No break here, as multiple heterocycles might be formed

                for rxn_type in COUPLING_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(rxn_type, rsmi):
                        coupling_reactions.append(
                            (depth, rxn_type, product_smiles, reactants_smiles)
                        )
                        if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break # Assuming only one coupling reaction type per step
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            traverse(child, new_depth)

    traverse(route)

    if not heterocycle_formations or not coupling_reactions:
        return False, findings_json
    else:
        # Structural constraint: both events must be present anywhere in the route
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "ring_formation",
                    "any_coupling_reaction"
                ],
                "context": "both events must be present anywhere in the route"
            }
        })

    earliest_heterocycle_depth = max([depth for depth, _ in heterocycle_formations])
    latest_coupling_depth = min([depth for depth, _, _, _ in coupling_reactions])

    timing_correct = earliest_heterocycle_depth > latest_coupling_depth

    if not timing_correct:
        return False, findings_json
    else:
        # Structural constraint: sequence
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "before": {
                    "target": "ring_formation",
                    "occurrence": "earliest"
                },
                "after": {
                    "target": "any_coupling_reaction",
                    "occurrence": "latest"
                }
            }
        })

    coupling_involves_heterocycle = False
    for _, _, product, reactants in coupling_reactions:
        product_has_heterocycle = False
        for heterocycle in final_product_heterocycles:
            if checker.check_ring(heterocycle, product):
                product_has_heterocycle = True
                break

        if product_has_heterocycle:
            reactant_has_heterocycle = False
            for reactant in reactants:
                for heterocycle in final_product_heterocycles:
                    if checker.check_ring(heterocycle, reactant):
                        reactant_has_heterocycle = True
                        break
                if reactant_has_heterocycle:
                    break

            if reactant_has_heterocycle:
                coupling_involves_heterocycle = True
                break

    if not coupling_involves_heterocycle:
        return False, findings_json
    else:
        # Structural constraint: coupling involves heterocycle as reactant
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "any_coupling_reaction",
                    "any_heterocycle_of_interest"
                ],
                "context": "the coupling reaction must have the heterocycle as a reactant"
            }
        })

    result = True
    return result, findings_json
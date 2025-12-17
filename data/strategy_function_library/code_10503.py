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
    # Standard cross-couplings
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
    "Negishi coupling",
    "Stille reaction_aryl",
    "Stille reaction_vinyl",
    "Stille reaction_benzyl",
    "Stille reaction_allyl",
    "Stille reaction_aryl OTf",
    "Stille reaction_vinyl OTf",
    "Stille reaction_benzyl OTf",
    "Stille reaction_allyl OTf",
    "Stille reaction_other",
    "Stille reaction_other OTf",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Aryllithium cross-coupling",
    "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl OTf",
    "Sonogashira acetylene_aryl OTf",
    "Sonogashira alkyne_alkenyl halide",
    "Sonogashira acetylene_alkenyl halide",
    "Sonogashira alkyne_alkenyl OTf",
    "Sonogashira acetylene_alkenyl OTf",
    "Sonogashira alkyne_acyl halide",
    "Sonogashira acetylene_acyl halide",
    "Buchwald-Hartwig",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Ullmann condensation",
    "decarboxylative_coupling",
    "Heck terminal vinyl",
    "Heck non-terminal vinyl",
    "Oxidative Heck reaction",
    "Oxidative Heck reaction with vinyl ester",
    # Amide couplings and acylations
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Identifies the late-stage incorporation of a nitrile-bearing fragment via a specific coupling reaction, where the nitrile group is preserved. This strategy highlights the chemoselective nature of the coupling, tolerating the nitrile functional group. The specific coupling reactions checked are defined in the `COUPLING_REACTIONS_OF_INTEREST` list.
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

    result = {"found": False, "depth": float("inf")}

    def dfs_traverse(node, depth=0):
        nonlocal result, findings_json
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a coupling reaction from the predefined list
                is_coupling = False
                found_coupling_rxn = None
                for rxn_type in COUPLING_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_coupling = True
                        found_coupling_rxn = rxn_type
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                        break

                # The strategy is only relevant for the specified coupling reactions.
                if is_coupling:
                    # Check for nitrile group in reactants
                    nitrile_reactants = []
                    nitrile_in_reactant_found = False
                    for reactant in reactants_smiles:
                        if checker.check_fg("Nitrile", reactant):
                            nitrile_reactants.append(reactant)
                            if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                            nitrile_in_reactant_found = True

                    if not nitrile_reactants:
                        return

                    # Verify nitrile is preserved in the product
                    nitrile_in_product_found = False
                    if checker.check_fg("Nitrile", product_smiles):
                        if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                        nitrile_in_product_found = True

                    if nitrile_in_reactant_found and nitrile_in_product_found:
                        # Update result if this is a shallower depth (later step) than previously found
                        if depth < result["depth"]:
                            result["found"] = True
                            result["depth"] = depth
                            # Add structural constraint if all conditions are met
                            findings_json["structural_constraints"].append({
                                "type": "co-occurrence",
                                "details": {
                                    "scope": "single_step",
                                    "description": "A reaction from a predefined list of coupling reactions must occur, and a Nitrile functional group must be present in both at least one reactant and the product of that same reaction step.",
                                    "targets": [
                                        "any_reaction_from_list",
                                        "Nitrile"
                                    ]
                                }
                            })

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)

    return result["found"], findings_json

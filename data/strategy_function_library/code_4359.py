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


CROSS_COUPLING_REACTIONS = [
    # Suzuki couplings
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with sulfonic esters",
    "{Suzuki}",
    # Stille couplings
    "Stille reaction_aryl",
    "Stille reaction_aryl OTf",
    "Stille reaction_vinyl",
    "Stille reaction_vinyl OTf",
    "Stille reaction_benzyl",
    "Stille reaction_benzyl OTf",
    "Stille reaction_allyl",
    "Stille reaction_allyl OTf",
    "Stille reaction_other",
    "Stille reaction_other OTf",
    "{Stille}",
    # Other cross-couplings
    "Negishi coupling",
    "Heck terminal vinyl",
    "Ullmann condensation",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
    "Aryllithium cross-coupling",
    "{Buchwald-Hartwig}",
    "{Negishi}",
    "{Heck_terminal_vinyl}",
    "{Heck_non-terminal_vinyl}",
    "decarboxylative_coupling",
]

AROMATIC_RINGS_OF_INTEREST = [
    "benzene", "pyridine", "furan", "thiophene", "pyrrole", "imidazole",
    "oxazole", "thiazole", "pyrazole", "naphthalene", "indole", "quinoline",
    "isoquinoline", "pyrimidine", "pyrazine", "triazole", "tetrazole",
    "benzothiophene", "benzoxazole", "benzimidazole", "benzotriazole",
    "indazole", "purine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis uses a late-stage aryl-aryl cross-coupling (like Suzuki coupling)
    as the final step to connect two aromatic/heteroaromatic fragments.
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

    # Track if we found a late-stage aryl-aryl coupling
    found_late_stage_coupling = False

    # Determine if the root node is a molecule or reaction
    root_is_molecule = route["type"] == "mol"

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_coupling, findings_json

        # For a molecule root, the first reaction is at depth 1
        # For a reaction root, the first reaction is at depth 0
        is_late_stage = (root_is_molecule and depth == 1) or (not root_is_molecule and depth == 0)

        if node["type"] == "reaction" and is_late_stage:
            print(f"Examining potential late-stage reaction at depth {depth}")
            # Add positional constraint if this is a late-stage reaction
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "cross_coupling",
                    "position": "last_stage"
                }
            })

            # Check if this is a cross-coupling reaction
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                parts = rsmi.split(">")
                reactants_str = parts[0]
                product_str = parts[-1]
                reactants = reactants_str.split(".")

                print(f"Reaction SMILES: {rsmi}")
                print(f"Number of reactants: {len(reactants)}")

                # Check for known cross-coupling reactions directly
                reaction_found = False
                for reaction_type in CROSS_COUPLING_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found late-stage {reaction_type}")
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        reaction_found = True
                        # No break here, as multiple reaction types might match, and we want to record all

                if reaction_found:
                    # Verify both reactants have aromatic rings
                    aromatic_reactants_count = 0
                    reactant_aromatic_rings_found = []
                    for r in reactants:
                        try:
                            reactant_has_aromatic_ring = False
                            for ring in AROMATIC_RINGS_OF_INTEREST:
                                if checker.check_ring(ring, r):
                                    reactant_aromatic_rings_found.append(ring)
                                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                                    aromatic_reactants_count += 1
                                    print(f"Reactant contains {ring} ring: {r}")
                                    reactant_has_aromatic_ring = True
                                    break # Found one aromatic ring in this reactant, move to next reactant
                        except Exception as e:
                            print(f"Error checking aromatic rings in reactant: {e}")

                    # Verify product has aromatic rings
                    has_product_aromatic = False
                    product_aromatic_rings_found = []
                    try:
                        for ring in AROMATIC_RINGS_OF_INTEREST:
                            if checker.check_ring(ring, product_str):
                                product_aromatic_rings_found.append(ring)
                                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                                has_product_aromatic = True
                                print(f"Product contains {ring} ring: {product_str}")
                                break
                    except Exception as e:
                        print(f"Error checking aromatic rings in product: {e}")

                    if aromatic_reactants_count >= 2:
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "aromatic_reactants",
                                "operator": ">=",
                                "value": 2
                            }
                        })

                    if aromatic_reactants_count >= 2 and has_product_aromatic:
                        print("Confirmed late-stage aryl-aryl coupling with aromatic fragments")
                        found_late_stage_coupling = True
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "multiple_aromatic_reactants",
                                    "aromatic_product"
                                ]
                            }
                        })
                        return # Found the desired pattern, no need to continue DFS

        # Determine the new depth for the recursive call
        new_depth = depth
        if node["type"] != "reaction": # If current node is 'chemical' (or 'mol'), depth increases when going to a reaction
            new_depth = depth + 1

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    return found_late_stage_coupling, findings_json

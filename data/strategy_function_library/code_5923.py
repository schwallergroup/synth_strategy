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


COUPLING_REACTIONS_FOR_CONVERGENCE = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Amide formation",
    "Esterification of Carboxylic Acids",
    "Williamson Ether Synthesis",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl halide",
    "Heck terminal vinyl",
    "Stille reaction_aryl",
    "Negishi coupling",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Ullmann condensation",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis strategy where two or more complex fragments (defined as having >10 atoms or at least one ring) are joined in the final step. The function specifically checks for common coupling reactions, including Suzuki, amide formation, and N-arylation, as defined in the COUPLING_REACTIONS_FOR_CONVERGENCE list.
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

    final_step_has_multiple_reactants = False
    debug_mode = True  # Set to True to enable debug prints

    def debug_print(message):
        if debug_mode:
            print(message)

    def dfs_traverse(node, depth=0):
        nonlocal final_step_has_multiple_reactants, findings_json

        debug_print(f"Traversing node at depth {depth}, type: {node['type']}")

        # Process reaction nodes - the final reaction step is at depth 1 in retrosynthetic traversal
        if node["type"] == "reaction" and depth == 1:  # Final reaction step
            debug_print(f"Examining final reaction step at depth {depth}")

            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                debug_print(f"Reaction SMILES: {rsmi}")

                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                debug_print(f"Found {len(reactants)} reactants")

                # Check if there are at least 2 reactants
                if len(reactants) >= 2:
                    # Record the structural constraint for multiple reactants
                    # This corresponds to the 'count' constraint in the strategy JSON
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "complex_reactants_in_final_step",
                            "operator": ">=",
                            "value": 2
                        }
                    })

                    # Check complexity of reactants
                    complex_reactants = []
                    for i, reactant in enumerate(reactants):
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                num_atoms = mol.GetNumAtoms()
                                num_rings = rdMolDescriptors.CalcNumRings(mol)
                                debug_print(
                                    f"Reactant {i+1}: {reactant} - Atoms: {num_atoms}, Rings: {num_rings}"
                                )

                                # Define complexity: has rings or is moderately large
                                if num_rings > 0 or num_atoms > 10:
                                    complex_reactants.append(reactant)
                        except Exception as e:
                            debug_print(f"Error processing reactant {i+1}: {e}")
                            continue

                    debug_print(f"Found {len(complex_reactants)} complex reactants")

                    # Check if there are at least 2 complex reactants
                    if len(complex_reactants) >= 2:
                        # Verify this is a coupling reaction
                        is_coupling = False
                        for rxn_type in COUPLING_REACTIONS_FOR_CONVERGENCE:
                            if checker.check_reaction(rxn_type, rsmi):
                                debug_print(f"Detected coupling reaction: {rxn_type}")
                                is_coupling = True
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                                break

                        # If it's a coupling reaction
                        if is_coupling:
                            final_step_has_multiple_reactants = True
                            # Record the structural constraint for positional (last_stage)
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": "convergent_coupling_reaction",
                                    "position": "last_stage"
                                }
                            })
                            debug_print(
                                "Detected convergent synthesis with complex fragments in final step"
                            )

        # Traverse children with incremented depth
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains the same when going from reaction to chemical
            new_depth = depth + 1 if node['type'] != 'reaction' else depth
            dfs_traverse(child, new_depth)

    # Start traversal from the root (final product)
    dfs_traverse(route)
    debug_print(f"Final result: {final_step_has_multiple_reactants}")

    return final_step_has_multiple_reactants, findings_json

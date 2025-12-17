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


SUZUKI_REACTION_NAMES = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with sulfonic esters",
    "Suzuki coupling with boronic esters OTf",
    "Suzuki coupling with boronic esters",
    "Suzuki",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a convergent synthesis strategy using late-stage Suzuki coupling
    to join two complex fragments.
    """
    print("Starting late_stage_suzuki_coupling_strategy analysis")
    # Track if we found a Suzuki coupling at depth 0 or 1 (final or penultimate step)
    found_late_stage_suzuki = False

    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_suzuki, findings_json

        if node["type"] == "reaction" and depth <= 1:  # Check depths 0 and 1
            print(f"Analyzing reaction at depth {depth}")
            # Check if this is a Suzuki coupling reaction
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    print(f"No reaction SMILES found in metadata at depth {depth}")
                    return

                print(f"Reaction SMILES at depth {depth}: {rsmi}")

                # Use the checker function to identify Suzuki coupling
                is_suzuki = False
                for reaction_name in SUZUKI_REACTION_NAMES:
                    if checker.check_reaction(reaction_name, rsmi):
                        print(f"Found {reaction_name} at depth {depth}")
                        is_suzuki = True
                        findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        break

                # If we found a Suzuki coupling, check if it's joining complex fragments
                if is_suzuki:
                    # Add positional constraint if Suzuki is found at depth <= 1
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Suzuki coupling",
                            "position": "late_stage (depth <= 1)"
                        }
                    })

                    # Get reactants and product
                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]
                    reactants = reactants_part.split(".")

                    # Check complexity of reactants (number of atoms as a simple metric)
                    complex_reactants = 0
                    for reactant in reactants:
                        if not reactant.strip():
                            continue

                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            atom_count = mol.GetNumAtoms()
                            print(f"Reactant has {atom_count} atoms: {reactant}")
                            if atom_count > 8:  # Consider fragments with >8 atoms as complex
                                complex_reactants += 1

                    # Verify product formation
                    product_mol = Chem.MolFromSmiles(product_part)
                    if not product_mol:
                        print(f"Could not parse product SMILES: {product_part}")

                    # If at least two complex reactants are being joined, it's a convergent synthesis
                    if complex_reactants >= 2 and product_mol:
                        print(
                            f"Found convergent synthesis with {complex_reactants} complex fragments at depth {depth}"
                        )
                        found_late_stage_suzuki = True
                        # Add count constraint if complex reactants condition is met
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "complex_reactants (>8 atoms) in Suzuki coupling",
                                "operator": ">=",
                                "value": 2
                            }
                        })
                    else:
                        print(
                            f"Suzuki coupling found at depth {depth} but not joining complex fragments or product invalid"
                        )

            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only from chemical to reaction node
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Completed analysis. Found late-stage Suzuki coupling: {found_late_stage_suzuki}")
    return found_late_stage_suzuki, findings_json

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
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Negishi coupling",
    "Stille reaction_aryl",
    "Heck terminal vinyl",
    "Sonogashira alkyne_aryl halide",
    "Buchwald-Hartwig",
    "Ullmann condensation",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis follows a convergent approach where two or more complex fragments
    are joined in a late-stage coupling reaction. It specifically checks for the following
    named reactions in the final synthetic step: Suzuki, Negishi, Stille, Heck, Sonogashira,
    Buchwald-Hartwig, and Ullmann condensation. A fragment is considered complex based on
    atom and ring counts.
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

    is_convergent = False

    def dfs_traverse(node, depth=0):
        nonlocal is_convergent, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # A late-stage reaction occurs at depth=1 (child of the root product at depth=0)
            if depth == 1:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if this is one of the specified coupling reactions
                is_coupling_found = False
                for name in COUPLING_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(name, rsmi):
                        is_coupling_found = True
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                        break

                if is_coupling_found:
                    # Record positional constraint
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "any_of_listed_coupling_reactions",
                            "position": "first_reaction_step"
                        }
                    })

                    # Count complex reactants (those with significant complexity)
                    complex_reactants = 0
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            atom_count = mol.GetNumAtoms()
                            ring_count = mol.GetRingInfo().NumRings()

                            # Define complexity based on multiple factors
                            is_complex = (atom_count > 12) and (ring_count >= 1 or atom_count > 15)

                            if is_complex:
                                complex_reactants += 1

                    # For convergent synthesis, we need at least 2 complex reactants
                    if complex_reactants >= 2:
                        is_convergent = True
                        # Record count constraint
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "complex_reactants",
                                "operator": ">=",
                                "value": 2
                            }
                        })

        # Continue traversing with incremented depth
        # Optimization to stop early is avoided to adhere to modification constraints
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth if current node is not a reaction
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return is_convergent, findings_json
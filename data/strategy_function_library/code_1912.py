from typing import Tuple, Dict, List
import copy
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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis uses a convergent approach, defined as a
    late-stage reaction (depth of 2 or less) where at least two complex
    fragments are joined. A fragment is considered complex if it has more than
    10 atoms or at least two rings.
    """
    findings_template = {
        "atomic_checks": {"named_reactions": [], "ring_systems": [], "functional_groups": []},
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    convergent_final_step = False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_final_step, findings_json

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            # The original code was getting depth from metadata, but the problem description
            # implies depth should be calculated during traversal. We'll use the passed depth.
            # depth = node.get("metadata", {}).get("depth", None)
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_part = rsmi.split(">")[0]

            # Check if this is a late-stage step (depth 0, 1, or 2)
            if depth is None or depth in [0, "0", 1, "1", 2, "2"]:
                # Record the positional constraint if met
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "convergent_reaction",
                        "position": "late_stage",
                        "max_depth": 2
                    }
                })
                reactants = reactants_part.split(".")

                # Count complex reactants (those with significant structure)
                complex_reactants = 0
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and (
                            mol.GetNumAtoms() > 10 or rdMolDescriptors.CalcNumRings(mol) >= 2
                        ):
                            complex_reactants += 1
                            print(
                                f"Complex reactant found: {reactant} with {mol.GetNumAtoms()} atoms and {rdMolDescriptors.CalcNumRings(mol)} rings"
                            )
                    except Exception as e:
                        print(f"Error processing reactant {reactant}: {e}")
                        continue

                # Determine if this is a convergent step
                if complex_reactants >= 2:
                    print(
                        f"Detected convergent synthesis with {complex_reactants} complex fragments at depth {depth}"
                    )
                    convergent_final_step = True
                    # Record the count constraint if met
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "complex_reactants_in_convergent_step",
                            "operator": ">=",
                            "value": 2
                        }
                    })

        # Continue traversing
        for child in node.get("children", []):
            # New depth calculation logic:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for its children (which are chemicals)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for its children (which are reactions)
                dfs_traverse(child, depth + 1)

    # Start traversal from the root with initial depth 0
    dfs_traverse(route, depth=0)

    if not convergent_final_step:
        print("No convergent synthesis pattern detected")

    return convergent_final_step, findings_json

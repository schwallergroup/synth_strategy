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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy involving late-stage N-alkylation
    to join complex fragments.
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

    has_late_stage_n_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_n_alkylation, findings_json

        # Focus on reactions at depths 0-1 (late stage)
        if node["type"] == "reaction" and depth <= 1:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if this is an N-alkylation reaction
                is_n_alkylation_primary = checker.check_reaction(
                    "N-alkylation of primary amines with alkyl halides", rsmi
                )
                is_n_alkylation_secondary = checker.check_reaction(
                    "N-alkylation of secondary amines with alkyl halides", rsmi
                )

                is_n_alkylation = is_n_alkylation_primary or is_n_alkylation_secondary

                if is_n_alkylation:
                    if depth <= 1:
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "N-alkylation",
                                "position": "late_stage"
                            }
                        })

                    if is_n_alkylation_primary:
                        findings_json["atomic_checks"]["named_reactions"].append("N-alkylation of primary amines with alkyl halides")
                    if is_n_alkylation_secondary:
                        findings_json["atomic_checks"]["named_reactions"].append("N-alkylation of secondary amines with alkyl halides")

                    # Check if we have at least 2 reactants (fragment joining)
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                    if len(reactant_mols) >= 2:
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "reactants_in_N-alkylation",
                                "operator": ">=",
                                "value": 2
                            }
                        })
                        # Count heavy atoms to ensure these are complex fragments
                        reactant_heavy_atoms = [
                            sum(1 for atom in r.GetAtoms() if atom.GetAtomicNum() > 1)
                            for r in reactant_mols
                            if r
                        ]
                        if any(count >= 10 for count in reactant_heavy_atoms):
                            print(
                                f"Detected late-stage N-alkylation joining complex fragments at depth {depth}"
                            )
                            has_late_stage_n_alkylation = True
                            findings_json["structural_constraints"].append({
                                "type": "count",
                                "details": {
                                    "target": "heavy_atoms_in_a_reactant",
                                    "operator": ">=",
                                    "value": 10
                                }
                            })

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when going from chemical to reaction.
            # Depth remains the same when going from reaction to chemical.
            new_depth = depth
            if node["type"] != "reaction":  # Current node is chemical
                new_depth = depth + 1
            # If current node is reaction, new_depth remains 'depth'

            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)
    return has_late_stage_n_alkylation, findings_json

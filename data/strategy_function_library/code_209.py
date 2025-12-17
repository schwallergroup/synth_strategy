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


PRESERVED_HETEROCYCLIC_SCAFFOLDS = [
    "benzothiophene",
    "benzoxazole",
    "benzimidazole",
    "indole",
    "benzofuran",
    "quinoline",
    "isoquinoline",
    "benzothiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy where a heterocyclic scaffold
    (like benzothiophene) is preserved throughout the synthesis while functional
    groups are modified.
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

    all_mols_have_scaffold = True
    mol_count = 0
    preserved_scaffold = None

    def dfs_traverse(node, depth=0):
        nonlocal all_mols_have_scaffold, mol_count, preserved_scaffold, findings_json

        if node["type"] == "mol":
            mol_count += 1
            mol_smiles = node["smiles"]

            mol = Chem.MolFromSmiles(mol_smiles)
            if not mol or mol.GetNumAtoms() <= 5:
                return

            has_scaffold = False
            for scaffold in PRESERVED_HETEROCYCLIC_SCAFFOLDS:
                if checker.check_ring(scaffold, mol_smiles):
                    has_scaffold = True
                    if scaffold not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(scaffold)
                    if preserved_scaffold is None:
                        preserved_scaffold = scaffold
                    elif preserved_scaffold != scaffold:
                        all_mols_have_scaffold = False
                        # Record inconsistency
                        findings_json["structural_constraints"].append({
                            "type": "negation",
                            "details": {
                                "target": "scaffold_type_is_inconsistent_across_route"
                            }
                        })
                    break

            if not has_scaffold:
                all_mols_have_scaffold = False
                # Record molecule without required scaffold
                findings_json["structural_constraints"].append({
                    "type": "negation",
                    "details": {
                        "target": "molecule_without_required_scaffold"
                    }
                })

        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # chemical node
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    strategy_valid = all_mols_have_scaffold and mol_count > 1 and preserved_scaffold is not None

    # Record structural constraints based on final flags
    if mol_count > 1:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "molecules",
                "operator": ">",
                "value": 1
            }
        })
    
    # The negation constraints are added during traversal if violated.
    # If they are not violated, they are implicitly met, but we only record violations.
    # The overall 'strategy_valid' captures the combination.

    return strategy_valid, findings_json

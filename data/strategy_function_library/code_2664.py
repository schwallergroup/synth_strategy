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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a diarylacetylene motif (an alkyne directly connecting two aromatic rings) is present in the final product and is carried over from at least one precursor molecule in the synthesis.
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

    # Helper function to check if a molecule has a diarylacetylene motif
    def has_aromatic_alkyne_connection(mol_smiles):
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return False

        # Direct connection: [aromatic_carbon]-[alkyne_carbon]#[alkyne_carbon]-[aromatic_carbon]
        pattern1 = Chem.MolFromSmarts("[c]-[#6]#[#6]-[c]")

        if mol.HasSubstructMatch(pattern1):
            # Add atomic checks for alkyne and aromatic ring if the motif is found
            if "alkyne" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("alkyne")
            if "diarylacetylene" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("diarylacetylene")
            if "aromatic ring" not in findings_json["atomic_checks"]["ring_systems"]:
                findings_json["atomic_checks"]["ring_systems"].append("aromatic ring")
            return True

        return False

    # The target molecule must have the motif for it to be "maintained"
    if route["type"] == "mol" and not has_aromatic_alkyne_connection(route["smiles"]):
        return False, findings_json

    # Track molecules with the alkyne linker at each depth
    molecules_with_alkyne = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            smiles = node["smiles"]

            if has_aromatic_alkyne_connection(smiles):
                if depth not in molecules_with_alkyne:
                    molecules_with_alkyne[depth] = []
                molecules_with_alkyne[depth].append(smiles)

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node['type'] == 'mol'
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if alkyne linker is found at multiple depths
    depths_with_alkyne = list(molecules_with_alkyne.keys())

    # For the alkyne to be maintained, it must be present in the target molecule (depth 0)
    # and at least one other precursor depth (> 0).
    if (
        depths_with_alkyne
        and min(depths_with_alkyne) == 0
        and max(depths_with_alkyne) > 0
    ):
        result = True
        # Add structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "diarylacetylene_in_final_product",
                    "diarylacetylene_in_precursor"
                ],
                "description": "The diarylacetylene motif must be present in the final product (depth 0) and also in at least one precursor molecule (depth > 0)."
            }
        })

    return result, findings_json

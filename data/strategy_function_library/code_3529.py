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
    Detects if a Boc-protected amine is present throughout the synthesis.
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

    # Track molecules with Boc groups at each depth
    molecules_with_boc = {}
    # Track all molecules at each depth
    all_molecules = {}
    # Track if we've seen a Boc deprotection reaction
    boc_deprotection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_deprotection_found, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Skip in-stock molecules (starting materials)
            if node.get("in_stock", False):
                return

            # Add this molecule to the appropriate depth
            if depth not in all_molecules:
                all_molecules[depth] = []
            all_molecules[depth].append(mol_smiles)

            # Check if this molecule has a Boc group
            if checker.check_fg("Boc", mol_smiles):
                if depth not in molecules_with_boc:
                    molecules_with_boc[depth] = []
                molecules_with_boc[depth].append(mol_smiles)
                if "Boc" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Boc")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Check if this is a Boc deprotection reaction
            if checker.check_reaction("Boc amine deprotection", rsmi):
                boc_deprotection_found = True
                if "Boc amine deprotection" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Boc amine deprotection")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            # Depth remains the same when traversing from reaction to chemical
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "mol"
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    result = True

    # If we found a Boc deprotection, Boc is not maintained throughout
    if boc_deprotection_found:
        result = False
        # Add the negation constraint if a deprotection was found
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "Boc amine deprotection",
                "scope": "route"
            }
        })

    # Check if all depths (except 0, which might be the final product) have a Boc group
    all_depths = set(all_molecules.keys())
    depths_with_boc = set(molecules_with_boc.keys())

    # Remove depth 0 if it exists (final product)
    if 0 in all_depths:
        all_depths.remove(0)

    # If we have no depths to check, return False
    if not all_depths:
        result = False # No relevant depths to check, so Boc is not 'throughout'

    # Check if Boc is present at all depths
    boc_throughout = all_depths.issubset(depths_with_boc)

    if result and boc_throughout: # Only add this constraint if the overall condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Boc",
                "position": "all_intermediate_stages"
            }
        })
    elif not boc_throughout and not boc_deprotection_found: # If Boc is not throughout and no deprotection was found, it means it was missing at some stage
        result = False

    return result, findings_json

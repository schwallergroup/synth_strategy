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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the N-dibenzyl group is maintained throughout
    the synthesis (no protection/deprotection of this group).
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

    # Track molecules in the main synthetic pathway
    main_pathway_mols = []

    def dfs_traverse(node, is_main_path=True):
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Only check molecules in the main synthetic pathway
            if is_main_path:
                main_pathway_mols.append(mol_smiles)
                print(f"Added to main pathway: {mol_smiles}")

            # If this is a leaf node (starting material), stop traversal
            if node.get("in_stock", False) or not node.get("children", []):
                return

        elif node["type"] == "reaction" and "children" in node:
            # For reaction nodes, the first child is typically the main product
            # in retrosynthetic direction
            for i, child in enumerate(node.get("children", [])):
                # First child is considered part of the main pathway
                # Others are reagents/catalysts
                dfs_traverse(child, is_main_path=(i == 0))
            return

        # Process children for molecule nodes
        for child in node.get("children", []):
            dfs_traverse(child, is_main_path)

    # Start traversal from root
    dfs_traverse(route)

    # Check if N-dibenzyl group is present in all main pathway molecules
    if not main_pathway_mols:
        print("No molecules found in main synthetic pathway")
        return False, findings_json

    # Check for N-dibenzyl in each molecule of the main pathway
    has_dibenzyl = []
    for i, mol_smiles in enumerate(main_pathway_mols):
        try:
            mol = Chem.MolFromSmiles(mol_smiles)
            if mol:
                # Check for benzyl groups attached to nitrogen
                dibenzyl_pattern = Chem.MolFromSmarts("[#7]([CH2]c1ccccc1)[CH2]c1ccccc1")
                matches = mol.GetSubstructMatches(dibenzyl_pattern)

                has_n_dibenzyl = len(matches) > 0
                has_dibenzyl.append(has_n_dibenzyl)

                if has_n_dibenzyl:
                    print(f"Found N-dibenzyl in molecule {i}: {mol_smiles}")
                    # Record finding: N-dibenzyl functional group
                    if "N-dibenzyl" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("N-dibenzyl")
                else:
                    print(f"No N-dibenzyl in molecule {i}: {mol_smiles}")
        except Exception as e:
            print(f"Error processing SMILES in dibenzyl detection: {e}")
            has_dibenzyl.append(False)

    # Check if all molecules in the main pathway have N-dibenzyl group
    all_have_dibenzyl = all(has_dibenzyl)

    if all_have_dibenzyl:
        print("N-dibenzyl group maintained throughout synthesis")
    else:
        missing_indices = [i for i, has_n_dibenzyl in enumerate(has_dibenzyl) if not has_n_dibenzyl]
        print(f"N-dibenzyl group missing in molecules at indices: {missing_indices}")
        # Record finding: negation (absence of N-dibenzyl on main pathway)
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "absence_of_N-dibenzyl_on_main_pathway"
            }
        })

    return all_have_dibenzyl, findings_json

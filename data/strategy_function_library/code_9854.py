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
    Detects if at least one Cbz-protected molecule is present at every depth of the synthesis. This is checked by identifying the Cbz group via a specific SMARTS pattern in all reactants and products.
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

    cbz_present_at_depths = set()
    max_depth = 0
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal cbz_present_at_depths, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            smiles = node.get("smiles", "")
            if smiles:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    benzyl_pattern = Chem.MolFromSmarts("[O;$(OC(=O)N)][CH2]c1ccccc1")
                    if mol and mol.HasSubstructMatch(benzyl_pattern):
                        cbz_present_at_depths.add(depth)
                        if "Carboxybenzyl group" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxybenzyl group")
                except Exception:
                    pass

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                benzyl_pattern = Chem.MolFromSmarts("[O;$(OC(=O)N)][CH2]c1ccccc1")

                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(benzyl_pattern):
                    cbz_present_at_depths.add(depth)
                    if "Carboxybenzyl group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Carboxybenzyl group")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(benzyl_pattern):
                        cbz_present_at_depths.add(depth + 1)
                        if "Carboxybenzyl group" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Carboxybenzyl group")
                        break
            except Exception:
                pass

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (chemicals)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (reactions)
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if Cbz is present at all depths
    all_depths = set(range(max_depth + 1))

    if all_depths and all_depths.issubset(cbz_present_at_depths):
        result = True
        findings_json["structural_constraints"].append({
            "type": "persistent_presence",
            "details": {
                "target": "Carboxybenzyl group",
                "condition": "Present at every depth of the synthesis."
            }
        })

    return result, findings_json
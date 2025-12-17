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


HETEROCYCLES_OF_INTEREST = [
    "furan",
    "pyrrole",
    "thiophene",
    "pyrazole",
    "imidazole",
    "oxazole",
    "thiazole",
    "triazole",
    "tetrazole",
    "pyridine",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "indole",
    "benzofuran",
    "benzothiophene",
    "benzimidazole",
    "benzoxazole",
    "benzothiazole",
    "quinoline",
    "isoquinoline",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects final products that contain a heterocycle-alkyne-aryl motif. The function specifically checks for the presence of a heterocycle from the HETEROCYCLES_OF_INTEREST list, a benzene or naphthalene ring, and an `aromatic-alkyne-aromatic` substructure.
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

    motif_present = False

    def find_final_product(node):
        # If it's a molecule and not in_stock
        if node["type"] == "mol" and not node.get("in_stock", False):
            # If it has no children or all children are reaction nodes
            if not node.get("children", []) or all(
                child["type"] == "reaction" for child in node.get("children", [])
            ):
                return node

        # Recursively check children
        for child in node.get("children", []):
            result = find_final_product(child)
            if result:
                return result

        return None

    def check_motif(mol_smiles):
        nonlocal findings_json
        mol = Chem.MolFromSmiles(mol_smiles)
        if not mol:
            return False

        # Check for heterocycle
        heterocycle_found = False
        for h in HETEROCYCLES_OF_INTEREST:
            if checker.check_ring(h, mol_smiles):
                findings_json["atomic_checks"]["ring_systems"].append(h)
                heterocycle_found = True
        if not heterocycle_found:
            return False

        # Check for aryl group
        has_aryl = False
        if checker.check_ring("benzene", mol_smiles):
            findings_json["atomic_checks"]["ring_systems"].append("benzene")
            has_aryl = True
        if checker.check_ring("naphthalene", mol_smiles):
            findings_json["atomic_checks"]["ring_systems"].append("naphthalene")
            has_aryl = True

        if not has_aryl:
            return False

        # Check connectivity: heterocycle-alkyne-aryl
        pattern = "[a;R]-!@[C]#[C]-!@[c;R]"  # Aromatic ring atom connected to alkyne connected to aromatic carbon
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))

        if bool(matches):
            findings_json["atomic_checks"]["functional_groups"].append("aromatic-alkyne-aromatic substructure")
        
        return bool(matches)

    # Find the final product
    final_product = find_final_product(route)
    if final_product:
        motif_present = check_motif(final_product['smiles'])
        if motif_present:
            # Add structural constraints if the motif is present
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "heterocycle-alkyne-aryl motif",
                    "position": "last_stage"
                }
            })
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "heterocycle_of_interest",
                        "aryl_of_interest",
                        "aromatic-alkyne-aromatic_substructure"
                    ]
                }
            })

    return motif_present, findings_json

def dfs_traverse(node, depth=0, visited=None):
    if visited is None:
        visited = set()

    node_id = id(node)
    if node_id in visited:
        return
    visited.add(node_id)

    node['depth'] = depth

    for child in node.get('children', []):
        # New depth calculation logic
        if node['type'] == 'reaction':
            # If current node is reaction, depth remains the same for children
            dfs_traverse(child, depth, visited)
        else:
            # If current node is not reaction (e.g., chemical), depth increases
            dfs_traverse(child, depth + 1, visited)

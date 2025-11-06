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


SNAR_REACTION_NAMES = [
    "nucl_sub_aromatic_ortho_nitro",
    "nucl_sub_aromatic_para_nitro",
    "heteroaromatic_nuc_sub",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy where a nitro group is carried through
    multiple steps to serve as an activating group for a key S-N-Ar reaction. It specifically
    checks for named reactions including 'nucl_sub_aromatic_ortho_nitro',
    'nucl_sub_aromatic_para_nitro', and 'heteroaromatic_nuc_sub'.
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

    # Track molecules with nitro groups and their atom indices
    nitro_molecules = []
    nitro_reactions = []
    snar_reactions = []

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check for nitro group in molecule
            if checker.check_fg("Nitro group", mol_smiles):
                findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
                # Get atom indices of nitro groups
                nitro_indices = checker.get_fg_atom_indices("Nitro group", mol_smiles)
                nitro_molecules.append((depth, mol_smiles, nitro_indices))
                print(f"Found molecule with nitro group at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                product = rsmi.split(">")[-1]

                # Check if product has nitro group
                if checker.check_fg("Nitro group", product):
                    findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
                    nitro_reactions.append((depth, rsmi))
                    print(f"Found reaction with nitro group in product at depth {depth}")

                    # Check if this is an SNAr reaction using named reactions
                    is_snar = False
                    for rxn in SNAR_REACTION_NAMES:
                        if checker.check_reaction(rxn, rsmi):
                            is_snar = True
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            break

                    if is_snar:
                        snar_reactions.append((depth, rsmi))
                        print(f"Found SNAr reaction at depth {depth}: {rsmi}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is 'mol', depth increases
            next_depth = depth + 1

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if nitro group persists through multiple steps
    # and is used in at least one SNAr reaction
    nitro_persistence = len(nitro_reactions) >= 2 and len(snar_reactions) >= 1

    # Additional check: verify that the nitro group is present in the final product
    # and was present in early stages (indicating persistence)
    if nitro_molecules:
        depths = [d for d, _, _ in nitro_molecules]
        if not depths:
            nitro_persistence = False
        else:
            nitro_persistence = nitro_persistence and (min(depths) == 0) and (max(depths) >= 2)

    # Record structural constraints if conditions are met
    if len(nitro_reactions) >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction_producing_nitro_group",
                "operator": ">=",
                "value": 2
            }
        })
    if len(snar_reactions) >= 1:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "SNAr_reaction",
                "operator": ">=",
                "value": 1
            }
        })
    if nitro_molecules and (min(depths) == 0):
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Nitro group",
                "position": "last_stage"
            }
        })
    if nitro_molecules and (max(depths) >= 2):
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "Nitro group",
                "position": "present_at_depth_gte_2"
            }
        })

    print(f"Nitro group persistence: {nitro_persistence}")
    print(f"Found {len(nitro_molecules)} molecules with nitro groups")
    print(f"Found {len(nitro_reactions)} reactions with nitro groups in products")
    print(f"Found {len(snar_reactions)} SNAr reactions")

    return nitro_persistence, findings_json

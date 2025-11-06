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
    This function detects a synthetic strategy involving late-stage ring fusion to form a tricyclic system.
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

    has_late_stage_ring_fusion = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_ring_fusion, findings_json

        # Update node depth
        node["depth"] = depth

        if node["type"] == "reaction" and depth <= 2:  # Late stage (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Create RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and reactant_mols:
                    # Count rings properly
                    reactant_ring_count = sum(
                        [mol.GetRingInfo().NumRings() for mol in reactant_mols]
                    )
                    product_ring_count = product_mol.GetRingInfo().NumRings()

                    # Check if this is a ring fusion reaction forming a tricyclic system
                    if product_ring_count > reactant_ring_count and product_ring_count >= 3:

                        # Count atoms that are in multiple rings to verify fusion
                        product_ring_info = product_mol.GetRingInfo()
                        atoms_in_multiple_rings = 0
                        for atom_idx in range(product_mol.GetNumAtoms()):
                            ring_count = product_ring_info.NumAtomRings(atom_idx)
                            if ring_count > 1:
                                atoms_in_multiple_rings += 1

                        has_fused_rings = atoms_in_multiple_rings >= 2

                        # Check if reactants already had fused rings
                        reactant_fused_atoms = 0
                        for mol in reactant_mols:
                            ring_info = mol.GetRingInfo()
                            for atom_idx in range(mol.GetNumAtoms()):
                                if ring_info.NumAtomRings(atom_idx) > 1:
                                    reactant_fused_atoms += 1

                        # Verify that a new fusion is happening in this reaction
                        new_fusion = (
                            has_fused_rings and atoms_in_multiple_rings > reactant_fused_atoms
                        )

                        if new_fusion:
                            has_late_stage_ring_fusion = True
                            # Record atomic check: named_reactions - ring_fusion
                            if "ring_fusion" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_fusion")
                            
                            # Record structural constraint: positional - late_stage
                            findings_json["structural_constraints"].append({
                                "type": "positional",
                                "details": {
                                    "target": "ring_fusion",
                                    "position": "late_stage",
                                    "max_depth": 2
                                }
                            })
                            # Record structural constraint: count - rings_in_product_of_ring_fusion
                            findings_json["structural_constraints"].append({
                                "type": "count",
                                "details": {
                                    "target": "rings_in_product_of_ring_fusion",
                                    "operator": ">=",
                                    "value": 3
                                }
                            })

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is chemical, depth increases
            next_depth = depth + 1

        # Traverse children with the determined next_depth
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)
    return has_late_stage_ring_fusion, findings_json

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
    Detects reactions that simultaneously form one or more rings and create at least one new stereocenter.
    This is identified by comparing the total ring count and chiral center count of the reactants versus the product.
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

    stereoselective_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal stereoselective_formation, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                    product = Chem.MolFromSmiles(product_smiles)

                    if product and all(r for r in reactants):
                        # Count stereocenters in reactants and product, including unassigned ones
                        reactant_chiral_centers = sum(
                            len(Chem.FindMolChiralCenters(r, includeUnassigned=True)) for r in reactants
                        )
                        product_chiral_centers = len(Chem.FindMolChiralCenters(product, includeUnassigned=True))

                        # Count rings in reactants and product
                        reactant_ring_count = sum(r.GetRingInfo().NumRings() for r in reactants)
                        product_ring_count = product.GetRingInfo().NumRings()

                        ring_formed = False
                        stereocenter_formed = False

                        if product_ring_count > reactant_ring_count:
                            ring_formed = True
                            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                        if product_chiral_centers > reactant_chiral_centers:
                            stereocenter_formed = True
                            if "stereocenter_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("stereocenter_formation")

                        # Check if new stereocenters were formed along with new rings
                        if ring_formed and stereocenter_formed:
                            stereoselective_formation = True
                            # Add the structural constraint if both conditions are met
                            structural_constraint_obj = {
                                "type": "co-occurrence",
                                "details": {
                                    "targets": [
                                        "ring_formation",
                                        "stereocenter_formation"
                                    ],
                                    "scope": "reaction"
                                }
                            }
                            if structural_constraint_obj not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(structural_constraint_obj)

                except Exception:
                    # Silently ignore errors in SMILES parsing
                    pass

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)
    return stereoselective_formation, findings_json
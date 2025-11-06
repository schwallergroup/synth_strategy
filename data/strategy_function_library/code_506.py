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
    This function detects if the synthetic route maintains a cyclopropyl group throughout the synthesis.
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

    # Track if cyclopropyl is present in the target molecule
    target_has_cyclopropyl = False

    # Track if the cyclopropyl group is maintained throughout the synthesis
    cyclopropyl_maintained = True

    # Track molecules that should have cyclopropyl based on synthetic path
    molecules_with_cyclopropyl = set()

    # Track reactions where cyclopropyl is created
    cyclopropyl_creation_reactions = set()

    def dfs_traverse(node, depth=0):
        nonlocal target_has_cyclopropyl, cyclopropyl_maintained, molecules_with_cyclopropyl, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_cyclopropyl = checker.check_ring("cyclopropane", mol_smiles)
            if has_cyclopropyl:
                if "cyclopropane" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("cyclopropane")

            # If this is the target molecule (depth 0), check if it has cyclopropyl
            if depth == 0:
                target_has_cyclopropyl = has_cyclopropyl

                # If target has cyclopropyl, add it to the set of molecules that should have it
                if has_cyclopropyl:
                    molecules_with_cyclopropyl.add(mol_smiles)
                else:
                    # If target doesn't have cyclopropyl, no need to check further
                    return

            # For non-target molecules, check if they should have cyclopropyl
            elif mol_smiles in molecules_with_cyclopropyl and not has_cyclopropyl:
                cyclopropyl_maintained = False

            # If this is a starting material, we don't need to check its children
            if node.get("in_stock", False):
                return

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # For reaction nodes, determine which reactants should have cyclopropyl
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            product = rsmi.split(">")[-1]
            reactants = rsmi.split(">")[0].split(".")

            # If the product should have cyclopropyl
            if product in molecules_with_cyclopropyl:
                # Check if this is a cyclopropyl-forming reaction
                product_has_cyclopropyl = checker.check_ring("cyclopropane", product)

                # Check if any reactant has cyclopropyl
                reactant_with_cyclopropyl = None
                for reactant in reactants:
                    if checker.check_ring("cyclopropane", reactant):
                        reactant_with_cyclopropyl = reactant
                        molecules_with_cyclopropyl.add(reactant)
                        break

                # If no reactant has cyclopropyl but product does, cyclopropyl was created in this reaction
                if not reactant_with_cyclopropyl and product_has_cyclopropyl:
                    cyclopropyl_creation_reactions.add(rsmi)
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                    # Since cyclopropyl was created here, previous molecules don't need to have it
                    # We don't add any reactants to molecules_with_cyclopropyl
                elif reactant_with_cyclopropyl:
                    # If a reactant has cyclopropyl, it should be maintained
                    pass

            # Check for ring destruction
            if not checker.check_ring("cyclopropane", product):
                for reactant in reactants:
                    if checker.check_ring("cyclopropane", reactant):
                        if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")
                        break

        # Traverse children
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "mol"
                dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Only return True if target has cyclopropyl AND cyclopropyl is maintained
    result = target_has_cyclopropyl and cyclopropyl_maintained

    if target_has_cyclopropyl:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "cyclopropane",
                "position": "last_stage"
            }
        })
    
    if not cyclopropyl_maintained:
        # This implies ring_destruction happened or cyclopropyl was lost
        # The negation constraint is met if cyclopropyl_maintained is False
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "ring_destruction"
            }
        })

    return result, findings_json

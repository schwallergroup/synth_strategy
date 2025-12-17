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


from rdkit import Chem

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a methoxy group is present in the final product and is preserved in all methoxy-containing precursors on the heuristically-determined main synthetic path.
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

    # Track if we've found a methoxy group in the final product
    methoxy_in_final_product = False
    # Track if the methoxy group is preserved throughout the synthesis
    methoxy_preserved = True

    def is_methoxy_present(smiles):
        """Check if a methoxy group is present in the molecule"""
        # Use the checker to specifically check for methoxy groups
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        # Check for methoxy group (OCH3)
        patt = Chem.MolFromSmarts("[O][C;H3]")
        if mol.HasSubstructMatch(patt):
            return True

        return False

    def dfs_traverse(node, depth=0, main_pathway_molecules=None):
        nonlocal methoxy_in_final_product, methoxy_preserved, findings_json

        if main_pathway_molecules is None:
            main_pathway_molecules = set()

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Check if this is the final product (depth 0)
            if depth == 0:
                methoxy_in_final_product = is_methoxy_present(mol_smiles)

                # If no methoxy in final product, strategy doesn't apply
                if not methoxy_in_final_product:
                    methoxy_preserved = False
                    return
                else:
                    # Record finding: methoxy in final product
                    if "methoxy" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("methoxy")
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "methoxy", "position": "last_stage"}})

                # Add final product to main pathway
                main_pathway_molecules.add(mol_smiles)

            # For intermediates (not starting materials), check if methoxy is preserved
            # Only check molecules that are part of the main synthetic pathway
            elif not node.get("in_stock", False) and mol_smiles in main_pathway_molecules:
                has_methoxy = is_methoxy_present(mol_smiles)

                if not has_methoxy and methoxy_in_final_product:
                    methoxy_preserved = False
                    # Record finding: loss of methoxy group
                    findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "loss_of_methoxy_group_on_main_pathway"}})

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            # For reaction nodes, identify the product and reactants
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # If the product is in the main pathway, add all reactants that have methoxy to the main pathway
            if product in main_pathway_molecules:
                for reactant in reactants:
                    if is_methoxy_present(reactant):
                        main_pathway_molecules.add(reactant)
                        if "methoxy" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("methoxy")

        # Determine the new depth for recursive calls based on the current node's type
        new_depth = depth
        if node["type"] != "reaction": # If current node is 'mol' (chemical), depth increases
            new_depth = depth + 1
        # If current node is 'reaction', depth remains the same

        # Traverse children with the determined new depth
        for child in node.get("children", []):
            dfs_traverse(child, new_depth, main_pathway_molecules)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is valid if methoxy is in final product and preserved throughout
    result = methoxy_in_final_product and methoxy_preserved
    return result, findings_json

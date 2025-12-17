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
    Detects if a bromine substituent is maintained throughout the entire synthesis.

    This function checks if bromine atoms are present in the main synthetic pathway,
    excluding starting materials that don't contribute bromine to the synthesis.
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

    # Track if we've found at least one molecule with bromine
    found_bromine = False
    # Track if bromine is maintained throughout synthesis
    maintains_bromine = True

    halide_fgs = [
        "Aromatic halide",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Alkenyl halide"
    ]

    def check_for_halide(smiles):
        has_halide = False
        for fg_name in halide_fgs:
            if checker.check_fg(fg_name, smiles):
                if fg_name not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append(fg_name)
                has_halide = True
        return has_halide

    def dfs_traverse(node, depth=0):
        nonlocal found_bromine, maintains_bromine

        # Check molecule nodes
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            has_bromine = check_for_halide(mol_smiles)

            # If this is the target molecule (depth 0), it must have bromine
            if depth == 0:
                if not has_bromine:
                    print(f"Target molecule doesn't have bromine: {mol_smiles}")
                    maintains_bromine = False
                else:
                    found_bromine = True
                    # Add structural constraint for target molecule having halide
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "any_halide_group",
                            "position": "last_stage",
                            "description": "The final product molecule must contain a halide functional group."
                        }
                    })

        # Check reaction nodes to ensure bromine transfers from reactants to products
        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant has bromine
            reactant_has_bromine = any(
                check_for_halide(r)
                for r in reactants
            )

            # Check if product has bromine
            product_has_bromine = check_for_halide(product)

            # If reactants have bromine but product doesn't, bromine was lost
            if reactant_has_bromine and not product_has_bromine:
                print(f"Bromine lost in reaction: {rsmi}")
                maintains_bromine = False
                # Add structural constraint for halide loss
                findings_json["structural_constraints"].append({
                    "type": "negation",
                    "details": {
                        "target": "halide_loss_in_reaction",
                        "description": "No reaction step is allowed to remove a halide functional group that was present in the reactants."
                    }
                })

            # If reactants have bromine, mark that we've found bromine in the synthesis
            if reactant_has_bromine:
                found_bromine = True

        # Continue traversing children
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

    # We need to have found at least one bromine and maintained it throughout
    result = found_bromine and maintains_bromine
    return result, findings_json

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
    Detects reactions where a molecule containing a cyclohexane scaffold is conserved from reactant to product, and the overall set of functional groups on the molecule changes.
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

    # Track if we found cyclohexane scaffold with modifications
    found_cyclohexane = False
    cyclohexane_modifications = 0

    def dfs_traverse(node, depth, max_depth):
        nonlocal found_cyclohexane, cyclohexane_modifications, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check if this molecule contains a cyclohexane ring
            if checker.check_ring("cyclohexane", mol_smiles):
                found_cyclohexane = True
                if "cyclohexane" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("cyclohexane")

        elif node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains cyclohexane
                if checker.check_ring("cyclohexane", product_smiles):
                    found_cyclohexane = True
                    if "cyclohexane" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("cyclohexane")

                    # Check for modifications on the cyclohexane
                    for reactant in reactants_smiles:
                        if checker.check_ring("cyclohexane", reactant):
                            # Both reactant and product have cyclohexane
                            # Canonicalize SMILES for comparison
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            product_mol = Chem.MolFromSmiles(product_smiles)

                            if reactant_mol and product_mol:
                                reactant_canonical = Chem.MolToSmiles(
                                    reactant_mol, isomericSmiles=True
                                )
                                product_canonical = Chem.MolToSmiles(
                                    product_mol, isomericSmiles=True
                                )

                                # Check if there's a structural change
                                if reactant_canonical != product_canonical:
                                    # Check for functional group changes as a proxy for modification
                                    product_fgs = []
                                    reactant_fgs = []

                                    fg_list = [
                                        "Primary alcohol", "Secondary alcohol", "Tertiary alcohol",
                                        "Ketone", "Aldehyde", "Ether", "Ester", "Carboxylic acid",
                                        "Primary amine", "Secondary amine", "Tertiary amine",
                                        "Alkene", "Alkyne", "Primary halide", "Secondary halide",
                                        "Tertiary halide",
                                    ]

                                    for fg in fg_list:
                                        if checker.check_fg(fg, product_smiles):
                                            product_fgs.append(fg)
                                            if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                                findings_json["atomic_checks"]["functional_groups"].append(fg)
                                        if checker.check_fg(fg, reactant):
                                            reactant_fgs.append(fg)
                                            if fg not in findings_json["atomic_checks"]["functional_groups"]:
                                                findings_json["atomic_checks"]["functional_groups"].append(fg)

                                    # If the set of functional groups differs, count it as a modification
                                    if set(product_fgs) != set(reactant_fgs):
                                        cyclohexane_modifications += 1
            except Exception as e:
                # Silently ignore errors in production
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth, max_depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1, max_depth)

    # This part is simplified for the example; a real implementation
    # would need to calculate max_depth for the route.
    # For now, we assume it's passed or calculated in the wrapper.
    max_depth_val = 10  # Placeholder
    dfs_traverse(route, 1, max_depth_val)

    # The strategy is present if we found cyclohexane with at least one modification
    strategy_present = found_cyclohexane and cyclohexane_modifications > 0

    if strategy_present:
        # Add the structural constraint if the strategy is present
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "cyclohexane_modification",
                "operator": ">",
                "value": 0
            }
        })

    return strategy_present, findings_json

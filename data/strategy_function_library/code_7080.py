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
    This function detects the strategy of late-stage introduction of fluorinated groups.
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

    # Track if we found fluorinated group addition in the last steps
    found_late_fluorination = False
    fluorination_steps = []

    # First pass to determine max depth for threshold calculation
    def get_max_depth(node, depth=0):
        max_d = depth
        for child in node.get("children", []):
            max_d = max(max_d, get_max_depth(child, depth + 1))
        return max_d

    max_depth = get_max_depth(route)
    # In retrosynthetic analysis, lower depths correspond to later stages in forward synthesis
    late_stage_threshold = max_depth * 0.6
    print(f"Max depth: {max_depth}, Late-stage threshold: {late_stage_threshold}")

    def count_fluorine_atoms(smiles):
        """Count the number of fluorine atoms in a molecule"""
        if not smiles:
            return 0
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return 0
        return len([atom for atom in mol.GetAtoms() if atom.GetSymbol() == "F"])

    def dfs_traverse(node, depth=0):
        nonlocal found_late_fluorination, fluorination_steps, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")
            product = product_part

            # Count fluorine atoms in reactants
            reactant_f_atoms = 0
            for reactant in reactants:
                if not reactant:
                    continue
                reactant_f_atoms += count_fluorine_atoms(reactant)

            # Count fluorine atoms in product
            product_f_atoms = count_fluorine_atoms(product)

            # Check if fluorinated groups were added
            fluorination_detected = False

            # Check by atom count
            if product_f_atoms > reactant_f_atoms:
                fluorination_detected = True
                print(f"Found fluorinated group addition at depth {depth}: +{product_f_atoms - reactant_f_atoms} F atoms")

            # Check for specific fluorination reactions
            fluorination_reaction_names = [
                "Aromatic fluorination",
                "Fluorination",
                "Hydrofluorination of methallyl alkenes",
                "Difluoroalkylation-Induced 1, 2-Heteroarene Migration of Allylic Alcohols"
            ]
            for r_name in fluorination_reaction_names:
                if checker.check_reaction(r_name, rsmi):
                    fluorination_detected = True
                    findings_json["atomic_checks"]["named_reactions"].append(r_name)
                    print(f"Confirmed fluorination reaction at depth {depth}")
                    break # Only add one instance of reaction if multiple match

            if fluorination_detected:
                fluorination_steps.append((depth, 1))
                # Check if this is a late-stage fluorination (lower depth in retrosynthesis)
                if depth <= late_stage_threshold:
                    found_late_fluorination = True
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "fluorination_event",
                            "position": "late_stage_retrosynthesis_le_60_percent_depth"
                        }
                    })
                    print(f"Late-stage fluorination detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage fluorination strategy detected: {found_late_fluorination}")
    return found_late_fluorination, findings_json

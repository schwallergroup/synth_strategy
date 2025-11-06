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


COUPLING_REACTIONS_OF_INTEREST = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Stille reaction_aryl",
    "Negishi coupling",
    "Heck terminal vinyl",
    "Sonogashira acetylene_aryl halide",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis follows a linear strategy without convergent steps.
    A linear strategy means each reaction step adds one new building block to the growing molecule,
    rather than combining multiple complex fragments.
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

    # Track if we found any convergent steps (more than one non-reagent reactant)
    found_convergent_step = False

    def is_likely_reagent(smiles):
        """
        Determine if a molecule is likely a reagent rather than a major reactant.
        Reagents include common solvents, bases, acids, catalysts, etc.
        """
        if not smiles:
            return True

        # Create RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Count heavy atoms (non-hydrogen)
        heavy_atom_count = mol.GetNumHeavyAtoms()

        # Small molecules are likely reagents
        if heavy_atom_count < 4:
            return True

        # Check for common reagent patterns
        if "[B]" in smiles and heavy_atom_count < 8:  # Boron reagents
            return True
        if "[Si]" in smiles and heavy_atom_count < 8:  # Silicon reagents
            return True
        if "P" in smiles and heavy_atom_count < 8:  # Phosphorus reagents
            return True

        return False

    def is_coupling_reaction(rsmi):
        """
        Check if the reaction is a coupling reaction, which is inherently convergent.
        """
        for reaction in COUPLING_REACTIONS_OF_INTEREST:
            if checker.check_reaction(reaction, rsmi):
                findings_json["atomic_checks"]["named_reactions"].append(reaction)
                return True

        return False

    def dfs_traverse(node, reaction, depth, max_depth):
        nonlocal found_convergent_step, findings_json

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Check if this is a coupling reaction
            if is_coupling_reaction(rsmi):
                found_convergent_step = True
                # Add the structural constraint for negation of named_reaction
                # This constraint is met if a coupling reaction IS found.
                # The overall strategy is linear if this constraint is NOT met.
                if {"type": "negation", "details": {"target": "named_reaction", "values": COUPLING_REACTIONS_OF_INTEREST, "description": "The route must not contain any of the specified coupling reactions, as they are inherently convergent."}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "named_reaction", "values": COUPLING_REACTIONS_OF_INTEREST, "description": "The route must not contain any of the specified coupling reactions, as they are inherently convergent."}})

            # Filter out molecules that are likely reagents
            non_reagent_reactants = [r for r in reactants_smiles if r and not is_likely_reagent(r)]

            if len(non_reagent_reactants) > 1:
                found_convergent_step = True
                # Add the structural constraint for negation of convergent_step
                # This constraint is met if a convergent step IS found.
                # The overall strategy is linear if this constraint is NOT met.
                if {"type": "negation", "details": {"target": "convergent_step", "description": "A reaction step must not have more than one non-reagent reactant."}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "convergent_step", "description": "A reaction step must not have more than one non-reagent reactant."}})

        # Traverse children
        for child in node.get("children", []):
            # The recursive call must be robust to missing reaction objects in non-reaction nodes
            child_reaction = child if child.get("type") == "reaction" else None
            
            # New depth calculation logic
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same

            dfs_traverse(child, child_reaction, new_depth, max_depth)

    # A full traversal is needed to determine max_depth, which is not performed here.
    # We pass placeholder values for depth and max_depth as they are not used by this function's logic.
    dfs_traverse(route, route, 0, 0)

    # Strategy is linear if no convergent steps were found
    is_linear = not found_convergent_step

    return is_linear, findings_json

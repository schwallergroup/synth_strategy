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


COUPLING_REACTIONS_FOR_CONVERGENCE = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Stille reaction_aryl",
    "Negishi coupling",
    "Heck terminal vinyl",
    "Sonogashira alkyne_aryl halide",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
    "Ullmann condensation",
    "Hiyama-Denmark Coupling",
    "Kumada cross-coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis strategy where two significant fragments
    are coupled in the middle stages of the synthesis.
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

    # Track if we found fragment coupling
    fragment_coupling_found = False
    # Track the depth at which fragment coupling occurs
    fragment_coupling_depth = None
    # Maximum depth in the route
    max_depth = 0

    def is_coupling_reaction(rsmi):
        # Check if the reaction is a known coupling reaction
        for rxn_type in COUPLING_REACTIONS_FOR_CONVERGENCE:
            if checker.check_reaction(rxn_type, rsmi):
                print(f"Detected coupling reaction: {rxn_type}")
                findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal fragment_coupling_found, fragment_coupling_depth, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Split reactants
                reactants = reactants_part.split(".")

                # Convert to molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product_part) if product_part else None

                if all(reactant_mols) and product_mol:
                    # Check if reactants are significant (at least 5 heavy atoms)
                    significant_reactants = [
                        mol for mol in reactant_mols if mol.GetNumHeavyAtoms() >= 5
                    ]

                    if len(significant_reactants) >= 2 and is_coupling_reaction(rsmi):
                        fragment_coupling_found = True
                        fragment_coupling_depth = depth
                        print(f"Fragment coupling detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if this is a mid-stage fragment coupling (in the middle half of synthesis)
    is_mid_stage = (
        fragment_coupling_depth is not None
        and fragment_coupling_depth >= (max_depth / 4)
        and fragment_coupling_depth <= (3 * max_depth / 4)
    )

    result = fragment_coupling_found and is_mid_stage
    if result:
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "fragment_coupling_reaction",
                "position_description": "The coupling reaction must occur in the middle half of the synthesis, at a depth (D) where max_depth/4 <= D <= 3*max_depth/4."
            }
        })

    print(f"Max depth: {max_depth}, Fragment coupling depth: {fragment_coupling_depth}")
    print(f"Convergent synthesis detection result: {result}")
    return result, findings_json
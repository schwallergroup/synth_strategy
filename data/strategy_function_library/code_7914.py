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


# Definition of the list of reactions for enumeration
FINAL_CYCLIZATION_REACTIONS = [
    "Paal-Knorr pyrrole synthesis",
    "Intramolecular transesterification/Lactone formation",
    "Intramolecular amination (heterocycle formation)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a linear fragment assembly strategy where fragments are added sequentially
    without convergent steps, building up to a final cyclization. The final cyclization
    is identified either by an increase in the number of rings or by checking for specific
    named reactions, including Paal-Knorr pyrrole synthesis, Intramolecular transesterification/Lactone formation, and Intramolecular amination (heterocycle formation).
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

    # Track the pattern
    is_linear = True
    has_final_cyclization = False
    step_count = 0
    first_step_processed = False

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, has_final_cyclization, step_count, first_step_processed, findings_json

        if node["type"] == "reaction":
            step_count += 1

            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if not reactant_mols or not product_mol:
                return

            # Check if this is a convergent step (more than 2 significant fragments combining)
            significant_fragments = 0
            for mol in reactant_mols:
                # Consider only non-trivial fragments (more than 5 atoms)
                if mol.GetNumAtoms() > 5:
                    significant_fragments += 1

            # In linear assembly, we should have at most 2 significant fragments
            if significant_fragments > 2:
                is_linear = False
                # Record the negation of convergent step if it occurs
                if {"type": "negation", "details": {"target": "convergent_step"}} not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "convergent_step"}})

            # Check for cyclization in final step (first reaction encountered in retrosynthetic traversal)
            if not first_step_processed:
                first_step_processed = True

                # Count rings in reactants and product
                reactant_ring_count = sum(len(Chem.GetSSSR(mol)) for mol in reactant_mols)
                product_ring_count = len(Chem.GetSSSR(product_mol))

                # Check if this is a cyclization reaction
                if product_ring_count > reactant_ring_count:
                    has_final_cyclization = True
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                elif any(
                    checker.check_reaction(rxn, rsmi)
                    for rxn in FINAL_CYCLIZATION_REACTIONS
                ):
                    has_final_cyclization = True
                    for rxn in FINAL_CYCLIZATION_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # From chemical to reaction
                new_depth = depth + 1
            # else: # From reaction to chemical, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Need at least 3 steps to be considered a meaningful linear strategy
    is_meaningful_linear = is_linear and step_count >= 3

    if step_count >= 3:
        if {"type": "count", "details": {"target": "reaction_step", "operator": ">=", "value": 3}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "reaction_step", "operator": ">=", "value": 3}})

    if is_linear:
        # If is_linear is True, it means no convergent step was found, so the negation constraint is met.
        # The negation constraint is added when is_linear becomes False, so if it's still True, we add it here.
        if {"type": "negation", "details": {"target": "convergent_step"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "convergent_step"}})

    if has_final_cyclization:
        if {"type": "positional", "details": {"target": "cyclization_event", "position": "last_stage"}} not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "cyclization_event", "position": "last_stage"}})

    result = is_meaningful_linear and has_final_cyclization

    return result, findings_json

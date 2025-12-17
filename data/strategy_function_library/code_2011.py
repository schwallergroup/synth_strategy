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


HETEROCYCLES_OF_INTEREST = [
    "benzoxazole",
    "benzothiazole",
    "benzimidazole",
    "quinoline",
    "isoquinoline",
    "indole",
    "purine",
    "pyrrole",
    "pyridine",
    "pyrazole",
    "imidazole",
    "oxazole",
    "thiazole",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "triazole",
    "tetrazole",
    "furan",
    "thiophene",
    "isoxazole",
    "isothiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """Detects a convergent synthesis where two distinct synthetic branches, each responsible for forming a specific heterocycle, are joined in a later-stage reaction. The list of heterocycles is defined in the HETEROCYCLES_OF_INTEREST constant."""
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    # Initialize tracking variables
    heterocycle_formations = []  # (depth, branch, ring, mol_smiles)
    convergent_reactions = []  # (depth, reaction_node, [branch_ids])

    # List of heterocyclic rings to check
    heterocycles = HETEROCYCLES_OF_INTEREST

    # Track molecule nodes to their branch IDs
    mol_to_branch = {}
    next_branch_id = 0

    def get_branch_id(mol_smiles):
        nonlocal next_branch_id, mol_to_branch
        if mol_smiles not in mol_to_branch:
            mol_to_branch[mol_smiles] = next_branch_id
            next_branch_id += 1
        return mol_to_branch[mol_smiles]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formations, convergent_reactions, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            current_branch = get_branch_id(mol_smiles)

            # Traverse children with the current branch ID
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)

        elif node["type"] == "reaction":
            # Extract reaction information
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Get branch IDs of reactants
            reactant_branches = []
            for r_smiles in reactants_smiles:
                if r_smiles in mol_to_branch:
                    reactant_branches.append(mol_to_branch[r_smiles])
                else:
                    # If reactant not seen before, assign a new branch ID
                    branch_id = get_branch_id(r_smiles)
                    reactant_branches.append(branch_id)

            # Check for heterocycle formation
            for ring in heterocycles:
                # Check if product contains heterocycle
                if checker.check_ring(ring, product_smiles):
                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                    # Check which reactants don't have the heterocycle
                    for i, r_smiles in enumerate(reactants_smiles):
                        if not checker.check_ring(ring, r_smiles):
                            branch_id = (
                                reactant_branches[i]
                                if i < len(reactant_branches)
                                else get_branch_id(r_smiles)
                            )
                            heterocycle_formations.append(
                                (depth, branch_id, ring, product_smiles)
                            )
                            # Mark ring formation as a named reaction
                            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

            # Check for convergent pattern - reaction with reactants from different branches
            unique_branches = set(reactant_branches)
            if len(unique_branches) >= 2:
                convergent_reactions.append((depth, node, list(unique_branches)))
                # Mark convergent reaction as a named reaction
                if "convergent_reaction" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("convergent_reaction")

            # Traverse children
            for child in node.get("children", []):
                # New logic: depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)

    # Start traversal from root
    dfs_traverse(route)

    # Get branches that form heterocycles
    heterocycle_branches = set(branch for _, branch, _, _ in heterocycle_formations)

    # Check if we have a convergent reaction that combines heterocycle-forming branches
    convergent_heterocycle_synthesis = False
    for depth, _, branches in convergent_reactions:
        # Check if this convergent reaction combines at least two heterocycle-forming branches
        heterocycle_branches_in_reaction = [b for b in branches if b in heterocycle_branches]
        if len(heterocycle_branches_in_reaction) >= 2:
            convergent_heterocycle_synthesis = True
            # Add structural constraints if the condition is met
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "heterocycle_forming_branch_in_convergent_reaction",
                    "operator": ">=",
                    "value": 2
                }
            })
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "ring_formation",
                    "after": "convergent_reaction"
                }
            })
            break

    return convergent_heterocycle_synthesis, findings_json
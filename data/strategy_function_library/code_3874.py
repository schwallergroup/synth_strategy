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


AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthesis route involves a late-stage (final step) amide bond formation
    between two complex fragments (>= 10 heavy atoms each), based on a curated list of named reactions.
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

    late_stage_amide_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_formation_found, findings_json

        # Process reaction nodes
        if node["type"] == "reaction":
            # Check if this is a late-stage reaction (depth <= 1)
            if depth <= 1 and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                # Skip if reaction is not bimolecular
                if len(reactants) < 2:
                    for child in node.get("children", []):
                        dfs_traverse(child, depth)
                    return

                # Check if this is an amide formation reaction
                is_amide_formation = False
                for reaction_type in AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_amide_formation = True
                        if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

                if is_amide_formation:
                    # Check if reactants are complex enough (at least 10 atoms each)
                    complex_reactants = 0
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.GetNumHeavyAtoms() >= 10:
                            complex_reactants += 1

                    # If at least 2 complex reactants, it's a late-stage amide formation
                    if complex_reactants >= 2:
                        late_stage_amide_formation_found = True
                        # Add structural constraints if found
                        if {"type": "positional", "details": {"target": "amide_formation", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "amide_formation", "position": "last_stage"}})
                        if {"type": "count", "details": {"target": "reactants_in_amide_formation_with_gte_10_heavy_atoms", "operator": ">=", "value": 2}} not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append({"type": "count", "details": {"target": "reactants_in_amide_formation_with_gte_10_heavy_atoms", "operator": ">=", "value": 2}})

        # Traverse children
        for child in node.get("children", []):
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return late_stage_amide_formation_found, findings_json
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


AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Carboxylic acid with primary amine to amide",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Schotten-Baumann_amide",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis with late-stage amide coupling.
    This means the final reaction combines two complex fragments via amide bond formation.
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

    is_convergent = False
    has_amide_coupling = False

    def is_complex_molecule(smiles):
        """Check if a molecule is complex enough to be considered a 'fragment'"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        # Consider a molecule complex if it has at least 10 atoms or contains a ring
        heavy_atom_count = mol.GetNumHeavyAtoms()
        has_ring = mol.GetRingInfo().NumRings() > 0

        return heavy_atom_count >= 10 or has_ring

    def check_reactant_complexity(reactants):
        """Check if at least two reactants are complex enough"""
        complex_reactants = [r for r in reactants if is_complex_molecule(r)]
        return len(complex_reactants) >= 2

    def dfs_traverse(node, depth=0):
        nonlocal is_convergent, has_amide_coupling, findings_json

        # Print current node for debugging
        if node["type"] == "mol":
            print(f"Depth {depth}, Molecule: {node['smiles'][:30]}...")
        else:
            print(f"Depth {depth}, Reaction node")

        # The final reaction is at depth 1 (one step back from target molecule)
        if node["type"] == "reaction" and depth == 1:
            # This is the 'late_stage' constraint
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "amide_coupling_reaction",
                    "position": "last_stage"
                }
            })

            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                print(f"Analyzing final reaction: {rsmi}")
                print(f"Reactants count: {len(reactants)}")

                # Check if there are at least 2 reactants and they're complex enough
                if len(reactants) >= 2:
                    if check_reactant_complexity(reactants):
                        is_convergent = True
                        # This is the 'count' constraint for complex reactants
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "complex_reactants_in_final_step",
                                "operator": ">=",
                                "value": 2
                            }
                        })
                    print(
                        f"Synthesis is convergent: {is_convergent} (2+ complex reactants in final step)"
                    )

                    # Check for amide coupling reactions
                    for reaction_name in AMIDE_COUPLING_REACTIONS:
                        if checker.check_reaction(reaction_name, rsmi):
                            has_amide_coupling = True
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                            print(f"Confirmed amide coupling reaction: {rsmi}")
                            break

                    print(f"Has amide coupling: {has_amide_coupling}")

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    result = is_convergent and has_amide_coupling

    if result:
        # This is the 'co-occurrence' constraint
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "late_stage_amide_coupling",
                    "convergent_final_step"
                ]
            }
        })

    print(
        f"Final result: {result} (Convergent: {is_convergent}, Amide coupling: {has_amide_coupling})
    )
    return result, findings_json

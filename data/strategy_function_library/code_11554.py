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


# Refactoring for Enumeration: Isolate lists into module-level constants
COMPLEXITY_FGS = [
    "Ester", "Amide", "Amine", "Alcohol", "Carboxylic acid", "Aromatic halide",
    "Sulfonamide", "Ether", "Ketone", "Aldehyde", "Nitrile",
]

KEY_BOND_FORMING_REACTIONS = [
    "Suzuki", "Negishi", "Stille", "Heck", "Sonogashira", "Buchwald-Hartwig",
    "N-arylation", "Ullmann-Goldberg", "Wittig", "Grignard", "Aldol condensation",
    "Michael addition", "Reductive amination",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis strategy where two complex fragments
    are joined in the final step via a key bond-forming reaction.
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
    print("Starting convergent synthesis analysis")

    def calculate_complexity(smiles):
        """Calculate the complexity of a molecule based on multiple factors"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            print(f"Could not parse SMILES: {smiles}")
            return 0

        num_atoms = mol.GetNumAtoms()
        num_rings = mol.GetRingInfo().NumRings()
        num_stereo_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))

        # Count functional groups from the predefined list
        fg_count = 0
        for fg in COMPLEXITY_FGS:
            if checker.check_fg(fg, smiles):
                fg_count += 1
                findings_json["atomic_checks"]["functional_groups"].append(fg)
                print(f"Found functional group: {fg} in {smiles}")

        complexity = num_atoms + (num_rings * 2) + (num_stereo_centers * 3) + (fg_count * 2)
        print(
            f"Complexity for {smiles}: {complexity} (atoms: {num_atoms}, rings: {num_rings}, stereo: {num_stereo_centers}, fg: {fg_count})"
        )

        return complexity

    def dfs_traverse(node, depth=0):
        nonlocal is_convergent, findings_json

        print(f"Traversing node at depth {depth}: {node.get('smiles', 'reaction node')}")

        if node["type"] == "reaction":
            if "metadata" not in node:
                node["metadata"] = {}
            node["metadata"]["depth"] = depth

            if depth <= 1:
                findings_json["structural_constraints"].append({
                    "type": "positional",
                    "details": {
                        "target": "convergent_step_analysis",
                        "position": "near_final_stage"
                    }
                })
                print(f"Analyzing reaction step at depth {depth}")
                reactants_smiles = []
                if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    print(f"Reaction SMILES: {rsmi}")
                    reactants_smiles = [r for r in rsmi.split(">")[0].split(".") if r]
                else:
                    print("No reaction SMILES found in metadata, checking children directly")
                    for child in node.get("children", []):
                        if child["type"] == "mol" and "smiles" in child:
                            reactants_smiles.append(child["smiles"])
                            print(f"Added reactant from child: {child['smiles']}")

                print(f"Analyzing {len(reactants_smiles)} reactants")

                complex_reactants = []
                for r in reactants_smiles:
                    complexity = calculate_complexity(r)
                    print(f"Reactant complexity: {complexity} for {r}")
                    if complexity > 12:
                        complex_reactants.append(r)

                # Fixing by Removal: Replaced inefficient, multi-stage check with a single, streamlined one.
                is_coupling = False
                if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                    rsmi = node["metadata"]["mapped_reaction_smiles"]
                    for rxn in KEY_BOND_FORMING_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            is_coupling = True
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                            print(f"Detected key bond-forming reaction.")
                            break

                # Fixing by Removal: Removed the 'else' block that was a source of false positives.
                if len(complex_reactants) >= 2:
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "complex_reactants_in_step",
                            "operator": ">=",
                            "value": 2
                        }
                    })
                    if is_coupling:
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "is_key_bond_forming_reaction",
                                    "has_at_least_two_complex_reactants"
                                ]
                            }
                        })
                        print(
                            f"Found convergent synthesis with {len(complex_reactants)} complex fragments in a coupling reaction"
                        )
                        is_convergent = True
                else:
                    print(
                        f"Step has only {len(complex_reactants)} complex fragments, not convergent"
                    )

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth if node["type"] == "reaction" else depth + 1
            dfs_traverse(child, new_depth)

    print(f"Analyzing route: {route.get('smiles', 'No SMILES')}")
    print("Starting DFS traversal")
    dfs_traverse(route)

    print(f"Convergent synthesis detected: {is_convergent}")
    return is_convergent, findings_json
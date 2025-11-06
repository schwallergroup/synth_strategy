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


COMPLEX_COUPLING_REACTION_TYPES = [
    "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters",
    "Suzuki coupling with boronic acids OTf",
    "Suzuki coupling with boronic esters OTf",
    "Negishi coupling",
    "Stille reaction_aryl",
    "Stille reaction_vinyl",
    "Heck terminal vinyl",
    "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl halide",
    "Ullmann condensation",
    "Buchwald-Hartwig",
    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
    "Kumada cross-coupling",
    "Hiyama-Denmark Coupling",
    "Ullmann-Goldberg Substitution amine",
    "Ullmann-Goldberg Substitution thiol",
    "Ullmann-Goldberg Substitution aryl alcohol",
    "Chan-Lam amine",
    "Chan-Lam alcohol",
    "Aryllithium cross-coupling",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthetic route involves a late-stage (depth <= 3) coupling of at least two complex fragments, where complexity is defined as a molecule containing two or more rings. The reaction is identified by matching against a predefined list of named cross-coupling reactions, such as Suzuki, Stille, Buchwald-Hartwig, and Negishi couplings.
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

    found_complex_coupling = False

    # Define the full strategy JSON for reference to pick structural constraints
    strategy_json = {
      "function_id": "code_11438",
      "filepath": "../data/merged_good_perf/code_11438.py",
      "description": "Detects if a synthetic route involves a late-stage (depth <= 3) coupling of at least two complex fragments, where complexity is defined as a molecule containing two or more rings. The reaction is identified by matching against a predefined list of named cross-coupling reactions, such as Suzuki, Stille, Buchwald-Hartwig, and Negishi couplings.",
      "atomic_checks": {
        "named_reactions": [
          "Suzuki coupling with boronic acids",
          "Suzuki coupling with boronic esters",
          "Suzuki coupling with boronic acids OTf",
          "Suzuki coupling with boronic esters OTf",
          "Negishi coupling",
          "Stille reaction_aryl",
          "Stille reaction_vinyl",
          "Heck terminal vinyl",
          "Sonogashira alkyne_aryl halide",
          "Sonogashira acetylene_aryl halide",
          "Ullmann condensation",
          "Buchwald-Hartwig",
          "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
          "Kumada cross-coupling",
          "Hiyama-Denmark Coupling",
          "Ullmann-Goldberg Substitution amine",
          "Ullmann-Goldberg Substitution thiol",
          "Ullmann-Goldberg Substitution aryl alcohol",
          "Chan-Lam amine",
          "Chan-Lam alcohol",
          "Aryllithium cross-coupling"
        ],
        "ring_systems": [],
        "functional_groups": []
      },
      "structural_constraints": [
        {
          "type": "count",
          "details": {
            "target": "reactants with 2 or more rings in a coupling reaction",
            "operator": ">=",
            "value": 2
          }
        },
        {
          "type": "positional",
          "details": {
            "target": "complex fragment coupling reaction",
            "position": "late_stage (depth <= 3)"
          }
        }
      ]
    }

    def dfs_traverse(node, depth=0):
        nonlocal found_complex_coupling, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            print(f"Checking reaction at depth {depth}: {rsmi}")

            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]
            reactants_smiles = reactants_part.split(".")

            # Check if this is a coupling reaction
            is_coupling_reaction = False
            for rxn_type in COMPLEX_COUPLING_REACTION_TYPES:
                if checker.check_reaction(rxn_type, rsmi):
                    is_coupling_reaction = True
                    if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    print(f"Identified as {rxn_type}")
                    break

            print(f"Is coupling reaction: {is_coupling_reaction}")

            if is_coupling_reaction and len(reactants_smiles) >= 2:
                # Convert to RDKit molecules
                reactant_mols = []
                for r in reactants_smiles:
                    mol = Chem.MolFromSmiles(r)
                    if mol is not None:
                        reactant_mols.append(mol)

                product_mol = Chem.MolFromSmiles(product_part)

                if len(reactant_mols) >= 2 and product_mol is not None:
                    try:
                        # Count rings in each reactant
                        ring_counts = []
                        for mol in reactant_mols:
                            ring_info = mol.GetRingInfo()
                            ring_count = ring_info.NumRings()
                            ring_counts.append(ring_count)

                        print(f"Ring counts in reactants: {ring_counts}")

                        # Check if at least two reactants have multiple rings
                        complex_fragments = sum(1 for count in ring_counts if count >= 2)

                        # Focus on reactions where at least two complex fragments are coupled
                        # Prioritize late-stage couplings (lower depth)
                        if complex_fragments >= 2 and depth <= 3:
                            print(f"Found complex fragment coupling at depth {depth}")
                            print(f"Reaction SMILES: {rsmi}")
                            print(f"Ring counts in reactants: {ring_counts}")
                            found_complex_coupling = True

                            # Add structural constraints to findings_json
                            # Constraint 1: reactants with 2 or more rings
                            constraint_1 = strategy_json["structural_constraints"][0]
                            if constraint_1 not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(copy.deepcopy(constraint_1))

                            # Constraint 2: late_stage (depth <= 3)
                            constraint_2 = strategy_json["structural_constraints"][1]
                            if constraint_2 not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(copy.deepcopy(constraint_2))

                    except Exception as e:
                        print(f"Error analyzing rings: {e}")

        # Continue traversing
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return found_complex_coupling, findings_json

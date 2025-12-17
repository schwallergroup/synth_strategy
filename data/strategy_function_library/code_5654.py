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


COMPLEX_FUSED_RINGS = [
    "naphthalene",
    "anthracene",
    "indole",
    "quinoline",
    "isoquinoline",
    "benzothiophene",
    "benzoxazole",
    "benzothiazole",
    "benzimidazole",
    "dibenzofuran",
    "dibenzothiophene",
    "carbazole",
    "acridine",
    "purine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the late-stage disconnection of a complex fused ring system, from a predefined list including naphthalene, indole, quinoline, etc.
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

    found_late_ring_disconnection = False
    
    # A robust pre-traversal to find the absolute max_depth is required.
    # This ensures the 'late-stage' definition is consistent.
    max_depth = 0
    q = [(route, 0)]
    visited = {id(route)}
    while q:
        node, depth = q.pop(0)
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            if id(child) not in visited:
                q.append((child, depth + 1))
                visited.add(id(child))

    def dfs_traverse(node, depth, max_depth):
        nonlocal found_late_ring_disconnection, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1].split(".")

                if not product_smiles or not reactants_smiles:
                    print("Skipping reaction with empty reactants or products")
                    return

                # Parse molecules
                reactant_mols = []
                for smi in reactants_smiles:
                    mol = Chem.MolFromSmiles(smi)
                    if mol:
                        reactant_mols.append(mol)
                    else:
                        print(f"Could not parse reactant SMILES: {smi}")

                product_mols = []
                for smi in product_smiles:
                    mol = Chem.MolFromSmiles(smi)
                    if mol:
                        product_mols.append(mol)
                    else:
                        print(f"Could not parse product SMILES: {smi}")

                if not reactant_mols or not product_mols:
                    return

                # Check if any reactant contains a complex ring system
                has_complex_ring_reactant = False
                for react_smi in reactants_smiles:
                    for ring in COMPLEX_FUSED_RINGS:
                        if checker.check_ring(ring, react_smi):
                            print(f"Found complex ring {ring} in reactant")
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)
                            has_complex_ring_reactant = True
                            break
                    if has_complex_ring_reactant:
                        break

                # Count rings in products and reactants
                product_rings = sum([len(Chem.GetSSSR(mol)) for mol in product_mols])
                reactant_rings = sum([len(Chem.GetSSSR(mol)) for mol in reactant_mols])

                print(f"Product rings: {product_rings}, Reactant rings: {reactant_rings}")

                # Check if a complex ring was disconnected in the late stage
                if reactant_rings > product_rings:
                    # This condition implies 'ring_destruction'
                    if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")

                    if has_complex_ring_reactant:
                        # This satisfies the 'co-occurrence' structural constraint
                        co_occurrence_constraint = {
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "ring_destruction",
                                    "complex_fused_ring_in_reactant"
                                ],
                                "definition": "A reaction step must both be a ring destruction (fewer rings in reactants than product) and have one of the complex fused rings present in a reactant molecule. The 'complex_fused_ring_in_reactant' target represents the check for any ring listed in atomic_checks.ring_systems."
                            }
                        }
                        if co_occurrence_constraint not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(co_occurrence_constraint)

                        # Define late stage as first third of synthesis depth
                        late_stage_threshold = max(1, max_depth // 3)
                        if depth <= late_stage_threshold:
                            print(
                                f"Detected complex ring disconnection at depth {depth}, threshold {late_stage_threshold}"
                            )
                            found_late_ring_disconnection = True
                            # This satisfies the 'positional' structural constraint
                            positional_constraint = {
                                "type": "positional",
                                "details": {
                                    "target": "disconnection_of_complex_fused_ring",
                                    "position": "late_stage_retrosynthesis",
                                    "definition": "The disconnection of a complex fused ring must occur within the first third of synthesis steps, starting from the final product (i.e., depth <= max_depth / 3)."
                                }
                            }
                            if positional_constraint not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(positional_constraint)

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth, max_depth)

    # Start traversal with pre-calculated max_depth
    dfs_traverse(route, 0, max_depth)

    print(f"Final result: {found_late_ring_disconnection}")
    return found_late_ring_disconnection, findings_json

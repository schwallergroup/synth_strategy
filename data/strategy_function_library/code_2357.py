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


from rdkit import Chem

# It is assumed the 'checker' module with 'check_ring' is available.
# import checker

COMMON_RINGS_OF_INTEREST = [
    "furan", "pyran", "pyrrole", "pyridine", "benzene", "naphthalene",
    "cyclopentane", "cyclohexane", "indole", "quinoline", "thiophene",
    "oxazole", "thiazole", "imidazole", "pyrazole", "triazole",
    "tetrazole", "piperidine", "piperazine", "morpholine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a late-stage cyclization strategy where a ring is formed
    in one of the final steps of the synthesis.
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

    cyclization_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal cyclization_detected, findings_json

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Focus on late-stage reactions (depth 0 or 1)
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    # Process each reactant separately
                    reactants_list = reactants_smiles.split(".")
                    reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_list]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if all(mol is not None for mol in reactants_mols) and product_mol is not None:
                        # Count rings in reactants and product
                        reactant_ring_count = sum(
                            len(mol.GetRingInfo().AtomRings()) for mol in reactants_mols
                        )
                        product_ring_count = len(product_mol.GetRingInfo().AtomRings())

                        # Check if a new ring was formed
                        if product_ring_count > reactant_ring_count:
                            # Verify it's a recognized ring structure
                            for ring_type in COMMON_RINGS_OF_INTEREST:
                                if checker.check_ring(ring_type, product_smiles) and not any(
                                    checker.check_ring(ring_type, r) for r in reactants_list
                                ):
                                    print(
                                        f"Late-stage cyclization detected at depth {depth}: {ring_type} ring formed"
                                    )
                                    cyclization_detected = True
                                    findings_json["atomic_checks"]["ring_systems"].append(ring_type)
                                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                                    
                                    # Structural constraints
                                    findings_json["structural_constraints"].append({
                                        "type": "positional",
                                        "details": {
                                            "target": "ring_formation",
                                            "position": "late_stage (depth <= 1)"
                                        }
                                    })
                                    findings_json["structural_constraints"].append({
                                        "type": "co-occurrence",
                                        "details": {
                                            "targets": [
                                                "ring_formation",
                                                "ring_system"
                                            ]
                                        }
                                    })
                                    findings_json["structural_constraints"].append({
                                        "type": "negation",
                                        "details": {
                                            "target": "presence of the newly formed ring system in reactants"
                                        }
                                    })
                                    break

                except Exception as e:
                    print(f"Error processing SMILES: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return cyclization_detected, findings_json

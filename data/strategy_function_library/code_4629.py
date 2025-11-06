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


HETEROCYCLIC_RINGS_OF_INTEREST = [
    "pyridine", "pyrrole", "furan", "thiophene", "imidazole", "oxazole",
    "thiazole", "pyrazole", "isoxazole", "isothiazole", "triazole",
    "tetrazole", "pyrimidine", "pyrazine", "pyridazine", "piperidine",
    "piperazine", "morpholine", "thiomorpholine", "oxirane", "aziridine",
    "oxetane", "azetidine", "tetrahydrofuran", "pyrrolidine",
    "tetrahydropyran", "oxepane", "azepane", "indole", "benzofuran",
    "benzothiophene", "benzimidazole", "benzoxazole", "benzothiazole",
    "quinoline", "isoquinoline",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the formation of a new heterocyclic ring from a predefined list (see HETEROCYCLIC_RINGS_OF_INTEREST)
    during the late stage (final or penultimate step) of the synthesis. The check confirms that a specific
    heterocycle type is present in the product but absent from all reactants.
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

    ring_formation_at_depth = None
    heterocyclic_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal ring_formation_at_depth, heterocyclic_formation, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if any heterocyclic ring is present in the product but not in reactants
                product_heterocycles = []
                for ring in HETEROCYCLIC_RINGS_OF_INTEREST:
                    if checker.check_ring(ring, product_smiles):
                        product_heterocycles.append(ring)

                reactant_heterocycles = []
                for r_smiles in reactants_smiles:
                    for ring in HETEROCYCLIC_RINGS_OF_INTEREST:
                        if checker.check_ring(ring, r_smiles):
                            reactant_heterocycles.append(ring)

                # Find heterocycles in product that aren't in reactants
                new_heterocycles = [
                    ring
                    for ring in product_heterocycles
                    if ring not in reactant_heterocycles
                ]

                if new_heterocycles:
                    heterocyclic_formation = True
                    ring_formation_at_depth = depth
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                    for ring_name in new_heterocycles:
                        if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                    print(
                        f"Heterocyclic ring formation detected at depth {depth}: {new_heterocycles}"
                    )

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (which are chemicals)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (which are reactions)
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    result = False
    # Check if heterocyclic ring formation occurred in late stage (depth 0 or 1)
    if (
        heterocyclic_formation
        and ring_formation_at_depth is not None
        and ring_formation_at_depth <= 1
    ):
        print("Late-stage heterocyclic ring formation strategy detected")
        result = True
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "ring_formation",
                "position": "late_stage"
            }
        })

    return result, findings_json

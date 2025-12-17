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


NITROGEN_HETEROCYCLES_OF_INTEREST = [
    "pyrrole", "pyridine", "pyrazole", "imidazole", "oxazole", "thiazole",
    "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole",
    "pyrrolidine", "piperidine", "piperazine", "morpholine", "indole",
    "quinoline", "isoquinoline", "benzimidazole", "benzoxazole", "benzothiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where a specific nitrogen-containing heterocycle is formed,
    followed by a subsequent reduction of a nitro group to an amine. The heterocycles
    checked are defined in the NITROGEN_HETEROCYCLES_OF_INTEREST list.
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

    # Track both conditions and their depths
    ring_formation_depth = float("inf")
    nitro_reduction_depth = float("inf")

    result = False

    def dfs_traverse(node, depth=0):
        nonlocal ring_formation_depth, nitro_reduction_depth, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for heterocyclic ring formation
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if all(reactant_mols) and product_mol:
                # Count rings in reactants and product
                reactant_ring_count = sum(mol.GetRingInfo().NumRings() for mol in reactant_mols)
                product_ring_count = product_mol.GetRingInfo().NumRings()

                # Check if a new ring is formed
                if product_ring_count > reactant_ring_count:
                    # Check if the product contains a nitrogen heterocycle
                    for heterocycle in NITROGEN_HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(heterocycle, product_smiles):
                            # Check if this heterocycle wasn't in any of the reactants
                            if not any(
                                checker.check_ring(heterocycle, r) for r in reactants_smiles
                            ):
                                print(
                                    f"Detected heterocyclic ring formation: {heterocycle} at depth {depth}"
                                )
                                ring_formation_depth = min(ring_formation_depth, depth)
                                if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                                if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                                break

            # Check for nitro reduction
            if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                print(f"Detected nitro reduction to amino group at depth {depth}")
                nitro_reduction_depth = min(nitro_reduction_depth, depth)
                if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # If current node is chemical, depth increases
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both conditions are met and in the correct sequence
    # Note: Lower depth means later stage in synthesis (closer to final product)
    if ring_formation_depth < float("inf") and nitro_reduction_depth < float("inf"):
        # Check if ring formation happens before nitro reduction in the synthetic route
        # (which means higher depth in the retrosynthetic tree)
        result = ring_formation_depth > nitro_reduction_depth
        if result:
            # Add structural constraints if both conditions are met and in correct sequence
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "ring_formation",
                        "Reduction of nitro groups to amines"
                    ]
                }
            })
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "Reduction of nitro groups to amines",
                    "after": "ring_formation"
                }
            })

    return result, findings_json

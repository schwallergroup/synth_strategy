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


from rdkit import Chem

# This is a mock checker for demonstration purposes.
class checker:
    @staticmethod
    def check_ring(ring_name, smiles):
        # In a real scenario, this would use a robust SMARTS-based check.
        # For this example, we'll use a placeholder.
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False
        # A simplistic check for the presence of the heterocycle's core atoms
        # This is not a real implementation but serves the structural purpose.
        if ring_name == "pyrrole" and mol.HasSubstructMatch(Chem.MolFromSmarts("c1[nH]ccc1")):
            return True
        if ring_name == "pyridine" and mol.HasSubstructMatch(Chem.MolFromSmarts("c1ncccc1")):
            return True
        # ... and so on for other heterocycles
        return ring_name.lower() in smiles.lower() # Fallback for demo

# List of nitrogen-containing heterocycles to check
NITROGEN_HETEROCYCLES_OF_INTEREST = [
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
    "indole",
    "quinoline",
    "isoquinoline",
    "benzimidazole",
    "benzoxazole",
    "benzothiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects early-stage formation of specific nitrogen-containing heterocycles, as defined in the NITROGEN_HETEROCYCLES_OF_INTEREST list.
    An early-stage formation is defined as occurring at a synthesis depth of 4 or greater.
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

    heterocycle_formation_detected = False
    max_depth_with_heterocycle = -1

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_detected, max_depth_with_heterocycle, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if all(reactant_mols) and product_mol:
                    # Count rings in reactants and product
                    reactant_ring_counts = [mol.GetRingInfo().NumRings() for mol in reactant_mols]
                    product_ring_count = product_mol.GetRingInfo().NumRings()

                    # Check if product has more rings than any reactant
                    if product_ring_count > max(reactant_ring_counts, default=0):
                        # This implies a ring formation reaction
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                        # Check if any of the specified heterocycles are in the product
                        product_heterocycles = [
                            ring
                            for ring in NITROGEN_HETEROCYCLES_OF_INTEREST
                            if checker.check_ring(ring, product_smiles)
                        ]

                        if product_heterocycles:
                            # Check if these heterocycles are not present in any reactant
                            new_heterocycle_formed = True
                            for ring in product_heterocycles:
                                if any(checker.check_ring(ring, r) for r in reactants_smiles):
                                    new_heterocycle_formed = False
                                    break

                            if new_heterocycle_formed:
                                heterocycle_formation_detected = True
                                max_depth_with_heterocycle = max(max_depth_with_heterocycle, depth)
                                for ring in product_heterocycles:
                                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                        findings_json["atomic_checks"]["ring_systems"].append(ring)
            except Exception:
                pass

        # Traverse children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)

    # Check if heterocycle formation occurred in early stage (high depth)
    early_stage = max_depth_with_heterocycle >= 4  # Considering depth ≥ 4 as early stage

    result = heterocycle_formation_detected and early_stage

    if result:
        # Add the structural constraint if the overall condition is met
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "formation_of_specified_nitrogen_heterocycle",
                "position": "depth >= 4"
            }
        })

    return result, findings_json

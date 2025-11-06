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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects formation of a fused thieno-pyridine heterocyclic system.
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

    thieno_pyridine_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal thieno_pyridine_formation_detected, findings_json

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product)

            has_thiophene_in_product = checker.check_ring("thiophene", product)
            has_pyridine_in_product = checker.check_ring("pyridine", product)

            if has_thiophene_in_product:
                if "thiophene" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("thiophene")
            if has_pyridine_in_product:
                if "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("pyridine")

            has_thieno_pyridine = has_thiophene_in_product and has_pyridine_in_product

            if has_thieno_pyridine and product_mol:
                thiophene_indices = checker.get_ring_atom_indices("thiophene", product)
                pyridine_indices = checker.get_ring_atom_indices("pyridine", product)

                if thiophene_indices and pyridine_indices:
                    fused_system = False
                    for thiophene_atoms in thiophene_indices:
                        for pyridine_atoms in pyridine_indices:
                            thiophene_set = set(thiophene_atoms)
                            pyridine_set = set(pyridine_atoms)

                            if thiophene_set.intersection(pyridine_set):
                                fused_system = True
                                break
                        if fused_system:
                            break

                    if fused_system:
                        # Record co-occurrence constraint
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    "thiophene",
                                    "pyridine"
                                ],
                                "scope": "product",
                                "note": "The strategy specifically checks for a fused thieno-pyridine system in the product."
                            }
                        })

                        has_fused_system_in_reactants = False
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                if checker.check_ring("thiophene", reactant) and checker.check_ring(
                                    "pyridine", reactant
                                ):
                                    r_thiophene_indices = checker.get_ring_atom_indices(
                                        "thiophene", reactant
                                    )
                                    r_pyridine_indices = checker.get_ring_atom_indices(
                                        "pyridine", reactant
                                    )

                                    if r_thiophene_indices and r_pyridine_indices:
                                        for r_thiophene_atoms in r_thiophene_indices:
                                            for r_pyridine_atoms in r_pyridine_indices:
                                                r_thiophene_set = set(r_thiophene_atoms)
                                                r_pyridine_set = set(r_pyridine_atoms)

                                                if r_thiophene_set.intersection(r_pyridine_set):
                                                    has_fused_system_in_reactants = True
                                                    break
                                            if has_fused_system_in_reactants:
                                                break
                            if has_fused_system_in_reactants:
                                break

                        if not has_fused_system_in_reactants:
                            thieno_pyridine_formation_detected = True
                            # Record negation constraint
                            findings_json["structural_constraints"].append({
                                "type": "negation",
                                "details": {
                                    "target": "fused_thieno-pyridine",
                                    "scope": "reactants",
                                    "note": "The strategy requires that the fused thieno-pyridine system is not present in the reactants of the formation step."
                                }
                            })
                            # Implicitly, if a new fused system is formed, it's a 'ring_formation' reaction
                            if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # If current node is reaction, depth remains same for children (chemical nodes)
                dfs_traverse(child, depth)
            else:
                # If current node is chemical, depth increases for children (reaction nodes)
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return thieno_pyridine_formation_detected, findings_json

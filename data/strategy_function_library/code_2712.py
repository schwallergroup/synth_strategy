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


HETEROCYCLES_OF_INTEREST = [
    "pyrrole",
    "pyridine",
    "pyrazole",
    "imidazole",
    "oxazole",
    "thiazole",
    "pyrimidine",
    "pyrazine",
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
    Checks for late-stage (depth <= 1) formation of specific heterocycles via a cyclization reaction that consumes a nitro group. The heterocycles of interest are defined in the HETEROCYCLES_OF_INTEREST list.
    """
    result = False
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    def dfs_traverse(node, depth=0):
        nonlocal result, findings_json

        if node["type"] == "reaction" and depth <= 1:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                reactant_list = reactants_smiles.split(".")
                reactants = [Chem.MolFromSmiles(r) for r in reactant_list if r]
                product = Chem.MolFromSmiles(product_smiles)

                if not all(reactants) or not product:
                    return

                has_nitro = False
                for r_smiles in reactant_list:
                    if checker.check_fg("Nitro group", r_smiles):
                        has_nitro = True
                        if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
                        break

                if has_nitro:
                    reactant_rings = sum(mol.GetRingInfo().NumRings() for mol in reactants)
                    product_rings = product.GetRingInfo().NumRings()

                    if product_rings > reactant_rings:
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                        heterocycle_found = False
                        for ring_name in HETEROCYCLES_OF_INTEREST:
                            if checker.check_ring(ring_name, product_smiles):
                                heterocycle_found = True
                                if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                                break

                        if heterocycle_found:
                            if not checker.check_fg("Nitro group", product_smiles):
                                if {"type": "negation", "details": {"target": "Nitro group in product"}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "negation", "details": {"target": "Nitro group in product"}})
                                
                                # All conditions met for the main result
                                result = True
                                
                                # Add structural constraints if not already present
                                if {"type": "positional", "details": {"target": "nitro_driven_cyclization", "position": "late_stage"}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "nitro_driven_cyclization", "position": "late_stage"}})
                                if {"type": "co-occurrence", "details": {"targets": ["Nitro group", "ring_formation", "specified_heterocycle"]}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["Nitro group", "ring_formation", "specified_heterocycle"]}})

            except Exception:
                pass

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return result, findings_json

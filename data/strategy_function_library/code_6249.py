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
    "isoxazole",
    "isothiazole",
    "oxadiazole",
    "thiadiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects C-N cross-coupling reactions (e.g., Buchwald-Hartwig, Ullmann) where one coupling partner is an amine and the other is a specific halogenated heterocycle. The heterocycle must be one of the following: pyridine, pyrazole, imidazole, oxazole, thiazole, pyrimidine, pyrazine, pyridazine, triazole, tetrazole, indole, quinoline, isoquinoline, benzimidazole, benzoxazole, benzothiazole, isoxazole, isothiazole, oxadiazole, or thiadiazole.
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

    found = False
    
    # Define the structural constraint object from the original strategy JSON
    structural_constraint_cn_cross_coupling = {
      "type": "co-occurrence",
      "details": {
        "description": "A C-N cross-coupling reaction must occur with a reactant that is both an aromatic halide and contains one of the specified heterocycles.",
        "targets": [
          "C-N cross-coupling (Buchwald-Hartwig type)",
          "Aromatic halide",
          "Specified N-heterocycle"
        ]
      }
    }

    def dfs_traverse(node, depth=0):
        nonlocal found, findings_json

        if node["type"] == "reaction" and not found:
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")

                try:
                    # Check if this is a Buchwald-Hartwig/Ullmann-Goldberg/N-arylation reaction
                    is_cn_cross_coupling = False
                    reaction_names = [
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                        "Buchwald-Hartwig"
                    ]
                    for r_name in reaction_names:
                        if checker.check_reaction(r_name, rsmi):
                            is_cn_cross_coupling = True
                            if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(r_name)
                            break

                    if is_cn_cross_coupling:
                        # Verify one of the reactants is a halogenated heterocycle from the list
                        has_halogen_heterocycle = False
                        for reactant in reactants:
                            if checker.check_fg("Aromatic halide", reactant):
                                if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                    findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")

                                for ring in HETEROCYCLES_OF_INTEREST:
                                    if checker.check_ring(ring, reactant):
                                        has_halogen_heterocycle = True
                                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                                        break
                            if has_halogen_heterocycle:
                                break

                        if has_halogen_heterocycle:
                            found = True
                            # Add the structural constraint if all conditions are met
                            if structural_constraint_cn_cross_coupling not in findings_json["structural_constraints"]:
                                findings_json["structural_constraints"].append(structural_constraint_cn_cross_coupling)

                except Exception:
                    pass

        for child in node.get("children", []):
            if not found:
                new_depth = depth
                if node["type"] != "reaction": # If current node is chemical, depth increases
                    new_depth += 1
                dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return found, findings_json

from typing import Tuple, Dict, List
import copy
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

AROMATIC_HETEROCYCLES_FOR_ALDEHYDE = [
    "pyridine",
    "pyrimidine",
    "pyrazine",
    "imidazole",
    "thiazole",
    "oxazole",
    "triazole",
    "furan",
    "thiophene",
    "isoxazole",
    "isothiazole",
]

SATURATED_N_HETEROCYCLES = [
    "piperidine",
    "pyrrolidine",
    "morpholine",
    "piperazine",
    "azetidine",
    "diazepane",
    "azepane",
    "pyrroline",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the use of three distinct building block types as reactants anywhere
    within the synthetic route: 1) an aldehyde-functionalized aromatic heterocycle
    from the list AROMATIC_HETEROCYCLES_FOR_ALDEHYDE, 2) a carbamate-protected
    saturated N-heterocycle from the list SATURATED_N_HETEROCYCLES, and 3) a
    benzyl halide.
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

    found_heterocycle = False
    found_protected_amine = False
    found_benzyl_fragment = False

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycle, found_protected_amine, found_benzyl_fragment, findings_json

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                for reactant_smiles in reactants_smiles:
                    if not reactant_smiles:
                        continue

                    # Check for heterocyclic aldehyde
                    if checker.check_fg("Aldehyde", reactant_smiles):
                        findings_json["atomic_checks"]["functional_groups"].append("Aldehyde")
                        for ring in AROMATIC_HETEROCYCLES_FOR_ALDEHYDE:
                            if checker.check_ring(ring, reactant_smiles):
                                findings_json["atomic_checks"]["ring_systems"].append(ring)
                                found_heterocycle = True
                                break # Found one relevant ring, no need to check others for this reactant

                    # Check for protected cyclic amine
                    if checker.check_fg("Carbamic ester", reactant_smiles):
                        findings_json["atomic_checks"]["functional_groups"].append("Carbamic ester")
                        for ring in SATURATED_N_HETEROCYCLES:
                            if checker.check_ring(ring, reactant_smiles):
                                findings_json["atomic_checks"]["ring_systems"].append(ring)
                                found_protected_amine = True
                                break # Found one relevant ring, no need to check others for this reactant

                    # Check for benzyl halide fragment
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol:
                        benzyl_pattern = Chem.MolFromSmarts("c1ccccc1C[!#1]")
                        if reactant_mol.HasSubstructMatch(benzyl_pattern):
                            findings_json["atomic_checks"]["ring_systems"].append("benzene") # Benzyl implies benzene ring
                            if (checker.check_fg("Primary halide", reactant_smiles)):
                                findings_json["atomic_checks"]["functional_groups"].append("Primary halide")
                                found_benzyl_fragment = True
                            if (checker.check_fg("Secondary halide", reactant_smiles)):
                                findings_json["atomic_checks"]["functional_groups"].append("Secondary halide")
                                found_benzyl_fragment = True

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is reaction, depth remains same for children (chemicals)
                dfs_traverse(child, depth)
            else:
                # If current node is chemical, depth increases for children (reactions)
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = found_heterocycle and found_protected_amine and found_benzyl_fragment

    if result:
        # Add the structural constraint if all three conditions are met
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "aldehyde-functionalized aromatic heterocycle",
                    "carbamate-protected saturated N-heterocycle",
                    "benzyl halide"
                ]
            }
        })

    # Deduplicate lists in findings_json
    for key in findings_json["atomic_checks"]:
        findings_json["atomic_checks"][key] = list(set(findings_json["atomic_checks"][key]))

    return result, findings_json

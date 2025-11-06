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


# Note: The checker and Chem objects are assumed to be globally available.
# from rdkit import Chem
# import checker

KNOWN_CYCLIZATION_REACTIONS = [
    "Paal-Knorr pyrrole synthesis",
    "Pictet-Spengler",
    "benzimidazole_derivatives_carboxylic-acid/ester",
    "benzimidazole_derivatives_aldehyde",
    "benzothiazole",
    "benzoxazole_arom-aldehyde",
    "benzoxazole_carboxylic-acid",
    "thiazole",
    "Niementowski_quinazoline",
    "tetrazole_terminal",
    "tetrazole_connect_regioisomere_1",
    "tetrazole_connect_regioisomere_2",
    "1,2,4-triazole_acetohydrazide",
    "1,2,4-triazole_carboxylic-acid/ester",
    "3-nitrile-pyridine",
    "spiro-chromanone",
    "pyrazole",
    "phthalazinone",
    "triaryl-imidazole",
    "Fischer indole",
    "Friedlaender chinoline",
    "benzofuran",
    "benzothiophene",
    "indole",
    "oxadiazole",
    "Formation of NOS Heterocycles",
    "Intramolecular amination of azidobiphenyls (heterocycle formation)",
    "Intramolecular amination (heterocycle formation)",
]

HETEROCYCLES_OF_INTEREST = [
    "furan", "pyran", "dioxane", "tetrahydrofuran", "tetrahydropyran",
    "oxirane", "oxetane", "oxolane", "oxane", "dioxolane", "pyrrole",
    "pyridine", "pyrazole", "imidazole", "oxazole", "thiazole",
    "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole",
    "pyrrolidine", "piperidine", "piperazine", "morpholine", "indole",
    "quinoline", "isoquinoline", "benzimidazole", "benzoxazole",
    "benzothiazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route involves a significant cyclization in the final step.
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

    late_cyclization_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal late_cyclization_detected, findings_json

        if node["type"] == "reaction" and depth <= 1:  # Final step
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                if product_smiles and all(r for r in reactants_smiles):
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if product_mol:
                        product_ring_info = product_mol.GetRingInfo()
                        product_ring_count = product_ring_info.NumRings()

                        reactant_ring_count = 0
                        for reactant in reactants_smiles:
                            if reactant.strip():
                                reactant_mol = Chem.MolFromSmiles(reactant)
                                if reactant_mol:
                                    reactant_ring_info = reactant_mol.GetRingInfo()
                                    reactant_ring_count += reactant_ring_info.NumRings()

                        if product_ring_count > reactant_ring_count:
                            is_cyclization_reaction = False
                            for rxn_type in KNOWN_CYCLIZATION_REACTIONS:
                                if checker.check_reaction(rxn_type, rsmi):
                                    is_cyclization_reaction = True
                                    if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                                    break

                            new_rings = []
                            for ring_name in HETEROCYCLES_OF_INTEREST:
                                if checker.check_ring(ring_name, product_smiles):
                                    ring_in_reactants = False
                                    for reactant in reactants_smiles:
                                        if checker.check_ring(ring_name, reactant):
                                            ring_in_reactants = True
                                            break
                                    if not ring_in_reactants:
                                        new_rings.append(ring_name)
                                        if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                                            findings_json["atomic_checks"]["ring_systems"].append(ring_name)

                            if is_cyclization_reaction or new_rings:
                                late_cyclization_detected = True
                                # Add structural constraint if detected
                                if {"type": "positional", "details": {"target": "ring_formation", "position": "last_stage"}} not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "ring_formation", "position": "last_stage"}})
            except Exception as e:
                pass

        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children (which are chemicals)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for children (which are reactions)
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_cyclization_detected, findings_json
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


HETEROCYCLIC_RINGS_OF_INTEREST = [
    "furan", "pyran", "dioxane", "tetrahydrofuran", "tetrahydropyran",
    "oxirane", "oxetane", "oxolane", "oxane", "dioxolane", "dioxolene",
    "trioxane", "dioxepane", "pyrrole", "pyridine", "pyrazole", "imidazole",
    "oxazole", "thiazole", "pyrimidine", "pyrazine", "pyridazine", "triazole",
    "tetrazole", "pyrrolidine", "piperidine", "piperazine", "morpholine",
    "thiomorpholine", "aziridine", "azetidine", "azepane", "diazepane",
    "indole", "quinoline", "isoquinoline", "purine", "carbazole", "acridine",
    "thiophene", "thiopyran", "thiirane", "thietane", "thiolane", "thiane",
    "dithiane", "dithiolane", "benzothiophene", "oxathiolane", "dioxathiolane",
    "thiazolidine", "oxazolidine", "isoxazole", "isothiazole", "oxadiazole",
    "thiadiazole", "benzoxazole", "benzothiazole", "benzimidazole", "pteridin",
    "phenothiazine", "phenoxazine", "dibenzofuran", "dibenzothiophene",
    "xanthene", "thioxanthene", "pyrroline", "pyrrolidone", "imidazolidine",
    "porphyrin", "indazole", "benzotriazole", "benzofuran",
]

HETEROCYCLE_FORMATION_REACTIONS_OF_INTEREST = [
    "benzofuran", "benzothiophene", "indole", "benzoxazole",
    "benzothiazole", "benzimidazole_derivatives_carboxylic-acid/ester",
    "benzimidazole_derivatives_aldehyde", "benzothiazole",
    "benzoxazole_arom-aldehyde", "benzoxazole_carboxylic-acid", "thiazole",
    "tetrazole_terminal", "tetrazole_connect_regioisomere_1",
    "tetrazole_connect_regioisomere_2", "1,2,4-triazole_acetohydrazide",
    "1,2,4-triazole_carboxylic-acid/ester", "3-nitrile-pyridine", "pyrazole",
    "oxadiazole", "imidazole", "Formation of NOS Heterocycles",
    "Paal-Knorr pyrrole synthesis", "Fischer indole", "Friedlaender chinoline",
    "Pictet-Spengler", "Niementowski_quinazoline",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition", "Pyrazole formation",
    "Azide-nitrile click cycloaddition to tetrazole",
    "Azide-nitrile click cycloaddition to triazole",
    "Intramolecular amination of azidobiphenyls (heterocycle formation)",
    "Intramolecular amination (heterocycle formation)",
    "Huisgen_Cu-catalyzed_1,4-subst",
    "Huisgen_Ru-catalyzed_1,5_subst",
    "Huisgen_disubst-alkyne",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects late-stage (depth <= 2) heterocycle formation. This is triggered if a reaction
    is a known heterocycle-forming named reaction (from the HETEROCYCLE_FORMATION_REACTIONS_OF_INTEREST list)
    OR if a heterocyclic ring (from the HETEROCYCLIC_RINGS_OF_INTEREST list) is present in the product
    but not in any of the reactants.
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

    found_late_stage_heterocycle = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_heterocycle, findings_json

        if node["type"] == "reaction" and depth <= 2:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                is_heterocycle_formation = False
                for reaction_type in HETEROCYCLE_FORMATION_REACTIONS_OF_INTEREST:
                    if checker.check_reaction(reaction_type, rsmi):
                        is_heterocycle_formation = True
                        if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                        break

                formed_heterocycles = []
                for ring in HETEROCYCLIC_RINGS_OF_INTEREST:
                    product_has_ring = checker.check_ring(ring, product_smiles)
                    if product_has_ring:
                        reactants_have_ring = any(
                            checker.check_ring(ring, r) for r in reactants_smiles if r
                        )
                        if not reactants_have_ring:
                            formed_heterocycles.append(ring)
                            if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring)

                if is_heterocycle_formation or formed_heterocycles:
                    found_late_stage_heterocycle = True
                    # Add the structural constraint if the condition is met
                    if {"type": "positional", "details": {"target": "heterocycle_formation", "position": "late_stage", "max_depth": 2}} not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "heterocycle_formation", "position": "late_stage", "max_depth": 2}})

            except Exception:
                pass

        if not found_late_stage_heterocycle:
            for child in node.get("children", []):
                # Determine the new depth based on the current node's type
                new_depth = depth
                if node["type"] != "reaction": # This means current node is 'chemical'
                    new_depth = depth + 1
                dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return found_late_stage_heterocycle, findings_json

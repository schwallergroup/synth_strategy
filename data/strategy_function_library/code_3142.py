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

HETEROCYCLE_TYPES = [
    "furan", "pyran", "dioxane", "tetrahydrofuran", "tetrahydropyran",
    "oxirane", "oxetane", "oxolane", "oxane", "dioxolane", "dioxolene",
    "pyrrole", "pyridine", "pyrazole", "imidazole", "oxazole", "thiazole",
    "pyrimidine", "pyrazine", "pyridazine", "triazole", "tetrazole",
    "pyrrolidine", "piperidine", "piperazine", "morpholine", "thiomorpholine",
    "indole", "quinoline", "isoquinoline", "benzoxazole", "benzothiazole",
    "benzimidazole", "purine", "carbazole", "acridine", "thiophene",
    "thiopyran", "thiirane", "thietane", "thiolane", "thiane", "dithiane",
    "dithiolane", "benzothiophene", "oxathiolane", "dioxathiolane",
    "thiazolidine", "oxazolidine", "isoxazole", "isothiazole",
    "oxadiazole", "thiadiazole",
]

HETEROCYCLE_FORMING_REACTIONS = [
    "Formation of NOS Heterocycles", "Paal-Knorr pyrrole synthesis",
    "benzothiazole formation from aldehyde", "benzothiazole formation from acyl halide",
    "benzothiazole formation from ester/carboxylic acid", "benzoxazole formation from aldehyde",
    "benzoxazole formation from acyl halide", "benzoxazole formation from ester/carboxylic acid",
    "benzoxazole formation (intramolecular)", "benzimidazole formation from aldehyde",
    "benzimidazole formation from acyl halide", "benzimidazole formation from ester/carboxylic acid",
    "tetrazole formation", "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition", "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Pyrazole formation", "Azide-nitrile click cycloaddition to tetrazole",
    "Azide-nitrile click cycloaddition to triazole",
    "Intramolecular amination of azidobiphenyls (heterocycle formation)",
    "Intramolecular amination (heterocycle formation)",
    "{benzimidazole_derivatives_carboxylic-acid/ester}", "{benzimidazole_derivatives_aldehyde}",
    "{benzothiazole}", "{benzoxazole_arom-aldehyde}", "{benzoxazole_carboxylic-acid}",
    "{thiazole}", "{Niementowski_quinazoline}", "{tetrazole_terminal}",
    "{tetrazole_connect_regioisomere_1}", "{tetrazole_connect_regioisomere_2}",
    "{Huisgen_Cu-catalyzed_1,4-subst}", "{Huisgen_Ru-catalyzed_1,5_subst}",
    "{Huisgen_disubst-alkyne}", "{1,2,4-triazole_acetohydrazide}",
    "{1,2,4-triazole_carboxylic-acid/ester}", "{3-nitrile-pyridine}", "{pyrazole}",
    "{Paal-Knorr pyrrole}", "{triaryl-imidazole}", "{Fischer indole}",
    "{Friedlaender chinoline}", "{benzofuran}", "{benzothiophene}", "{indole}",
    "{oxadiazole}", "{imidazole}",
    "1,2,4-oxadiazol-5(2H)-one synthesis from nitrile, hydrogen carbonate, and hydroxylamine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route involves formation of a heterocycle
    in the final step (late-stage heterocycle formation).
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

    final_heterocycle_formation = False

    def find_final_reaction(node, depth=0, path=None):
        """Find the final reaction in the synthetic route (closest to target)"""
        if path is None:
            path = []

        if depth == 0 and node["type"] == "mol":
            for child in node.get("children", []):
                if child["type"] == "reaction":
                    return child

        if node["type"] == "reaction":
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    has_reaction_children = False
                    for grandchild in child.get("children", []):
                        if grandchild["type"] == "reaction":
                            has_reaction_children = True
                            break

                    if not has_reaction_children:
                        return node

        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            result = find_final_reaction(child, new_depth, path + [node])
            if result:
                return result

        return None

    final_reaction = find_final_reaction(route)

    if final_reaction and "metadata" in final_reaction and "rsmi" in final_reaction["metadata"]:
        rsmi = final_reaction["metadata"]["rsmi"]

        try:
            parts = rsmi.split(">")
            if len(parts) >= 3:
                reactants_smiles = parts[0].split(".")
                product_smiles = parts[2]

                reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if product_mol and all(mol is not None for mol in reactants_mols):
                    for rxn_type in HETEROCYCLE_FORMING_REACTIONS:
                        if checker.check_reaction(rxn_type, rsmi):
                            final_heterocycle_formation = True
                            if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                            break

                    if not final_heterocycle_formation:
                        heterocycles_in_product = []
                        for heterocycle in HETEROCYCLE_TYPES:
                            if checker.check_ring(heterocycle, product_smiles):
                                heterocycles_in_product.append(heterocycle)
                                if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(heterocycle)

                        for heterocycle in heterocycles_in_product:
                            heterocycle_in_reactants = any(
                                checker.check_ring(heterocycle, r) for r in reactants_smiles
                            )

                            if not heterocycle_in_reactants:
                                final_heterocycle_formation = True
                                # The heterocycle was already added to findings_json["atomic_checks"]["ring_systems"]
                                # when it was found in the product.
                                break
        except Exception as e:
            pass

    if final_heterocycle_formation:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "heterocycle_formation",
                "position": "last_stage"
            }
        })

    return final_heterocycle_formation, findings_json

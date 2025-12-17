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


CYCLIZATION_REACTION_NAMES = [
    "Formation of NOS Heterocycles", "Paal-Knorr pyrrole synthesis",
    "Intramolecular transesterification/Lactone formation", "Diels-Alder",
    "Pauson-Khand reaction", "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition", "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Pyrazole formation", "Azide-nitrile click cycloaddition to tetrazole",
    "Azide-nitrile click cycloaddition to triazole", "Michael-induced ring closure from hydrazone",
    "Michael-induced ring closure from diazoalkane", "[3+2]-cycloaddition of hydrazone and alkyne",
    "[3+2]-cycloaddition of hydrazone and alkene", "[3+2]-cycloaddition of diazoalkane and alkyne",
    "[3+2]-cycloaddition of diazoalkane and alkene", "[3+2]-cycloaddition of diazoalkane and alpha-alkyne",
    "[3+2]-cycloaddition of diazoalkane and alpha-alkene",
    "Intramolecular amination of azidobiphenyls (heterocycle formation)",
    "Intramolecular amination (heterocycle formation)", "Pictet-Spengler",
    "benzimidazole_derivatives_carboxylic-acid/ester", "benzimidazole_derivatives_aldehyde",
    "benzothiazole", "benzoxazole_arom-aldehyde", "benzoxazole_carboxylic-acid",
    "thiazole", "Niementowski_quinazoline", "tetrazole_terminal",
    "tetrazole_connect_regioisomere_1", "tetrazole_connect_regioisomere_2",
    "Huisgen_Cu-catalyzed_1,4-subst", "Huisgen_Ru-catalyzed_1,5_subst",
    "Huisgen_disubst-alkyne", "1,2,4-triazole_acetohydrazide",
    "1,2,4-triazole_carboxylic-acid/ester", "3-nitrile-pyridine", "spiro-chromanone",
    "pyrazole", "phthalazinone", "Paal-Knorr pyrrole", "triaryl-imidazole",
    "Fischer indole", "Friedlaender chinoline", "benzofuran", "benzothiophene",
    "indole", "oxadiazole", "piperidine_indole", "imidazole",
]

ADDITIONAL_RING_FORMING_REACTIONS = [
    "Williamson Ether Synthesis (intra to epoxy)", "Diels-Alder (ON bond)",
    "Retro-Diels-Alder from oxazole", "Production of 2H,1-benzopyrans",
    "Benzothiazole formation from aldehyde", "Benzothiazole formation from acyl halide",
    "Benzothiazole formation from ester/carboxylic acid", "Benzoxazole formation from aldehyde",
    "Benzoxazole formation from acyl halide", "Benzoxazole formation from ester/carboxylic acid",
    "Benzoxazole formation (intramolecular)", "Benzimidazole formation from aldehyde",
    "Benzimidazole formation from acyl halide", "Benzimidazole formation from ester/carboxylic acid",
]

ALL_CYCLIZATION_REACTIONS = CYCLIZATION_REACTION_NAMES + ADDITIONAL_RING_FORMING_REACTIONS

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage cyclization strategy. A reaction is identified as a cyclization if it either results in a net increase in the number of rings or if its name matches a comprehensive list of known cyclization reactions (see ALL_CYCLIZATION_REACTIONS). The strategy is flagged if such a reaction occurs in the latter half of the synthesis (i.e., closer to the final product).
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

    cyclization_depth = None
    max_depth = 0
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal cyclization_depth, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                is_cyclization_rxn = False
                cyclization_rxn_name = None
                for rxn in ALL_CYCLIZATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        is_cyclization_rxn = True
                        cyclization_rxn_name = rxn
                        break

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if all(mol is not None for mol in reactant_mols) and product_mol is not None:
                    total_reactant_rings = sum(
                        mol.GetRingInfo().NumRings() for mol in reactant_mols
                    )
                    product_rings = product_mol.GetRingInfo().NumRings()

                    ring_increase = product_rings > total_reactant_rings

                    if ring_increase or is_cyclization_rxn:
                        if cyclization_depth is None or depth < cyclization_depth:
                            cyclization_depth = depth

                        # Record atomic check for named reaction if applicable
                        if is_cyclization_rxn and cyclization_rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(cyclization_rxn_name)
                        
                        # Record atomic check for ring formation if applicable
                        if ring_increase and "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children
                dfs_traverse(child, depth)
            else:
                # If current node is not a reaction (e.g., chemical), depth increases
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if cyclization_depth is not None and max_depth > 0:
        is_late_stage = cyclization_depth <= max_depth / 2
        if is_late_stage:
            result = True
            # Add structural constraints if the strategy is found
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "cyclization_event",
                    "operator": ">=",
                    "value": 1
                }
            })
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "cyclization_event",
                    "position": "latter_half (depth <= max_depth / 2)"
                }
            })

    return result, findings_json

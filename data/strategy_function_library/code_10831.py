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


HETEROCYCLIC_RINGS = [
    "furan",
    "pyran",
    "dioxane",
    "tetrahydrofuran",
    "tetrahydropyran",
    "oxirane",
    "oxetane",
    "oxolane",
    "oxane",
    "dioxolane",
    "dioxolene",
    "trioxane",
    "dioxepane",
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
    "pyrrolidine",
    "piperidine",
    "piperazine",
    "morpholine",
    "thiomorpholine",
    "aziridine",
    "azetidine",
    "azepane",
    "diazepane",
    "indole",
    "quinoline",
    "isoquinoline",
    "purine",
    "carbazole",
    "acridine",
    "thiophene",
    "thiopyran",
    "thiirane",
    "thietane",
    "thiolane",
    "thiane",
    "dithiane",
    "dithiolane",
    "benzothiophene",
    "oxathiolane",
    "dioxathiolane",
    "thiazolidine",
    "oxazolidine",
    "isoxazole",
    "isothiazole",
    "oxadiazole",
    "thiadiazole",
    "benzoxazole",
    "benzothiazole",
    "benzimidazole",
    "pteridin",
    "phenothiazine",
    "phenoxazole",
    "dibenzofuran",
    "dibenzothiophene",
    "xanthene",
    "thioxanthene",
    "pyrroline",
    "pyrrolidone",
    "imidazolidine",
    "porphyrin",
    "indazole",
    "benzotriazole",
]

HETEROCYCLE_FORMING_REACTIONS = [
    "Formation of NOS Heterocycles",
    "Paal-Knorr pyrrole synthesis",
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
    "Huisgen_Cu-catalyzed_1,4-subst",
    "Huisgen_Ru-catalyzed_1,5_subst",
    "1,2,4-triazole_acetohydrazide",
    "1,2,4-triazole_carboxylic-acid/ester",
    "3-nitrile-pyridine",
    "pyrazole",
    "Fischer indole",
    "Friedlaender chinoline",
    "benzofuran",
    "benzothiophene",
    "indole",
    "oxadiazole",
    "imidazole",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Pyrazole formation",
    "Azide-nitrile click cycloaddition to tetrazole",
    "Azide-nitrile click cycloaddition to triazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis strategy where at least two starting materials containing specific heterocycles are combined in at least one convergent reaction. The reaction is considered part of the strategy if it involves a heterocycle-containing reactant and either produces a product with a heterocycle or is a known heterocycle-forming reaction. The specific heterocycles and reaction types are defined in the module-level constants `HETEROCYCLIC_RINGS` and `HETEROCYCLE_FORMING_REACTIONS`.
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

    # Track heterocycle-containing starting materials
    heterocycle_starting_materials = []
    # Track convergent reactions that combine heterocycles
    convergent_heterocycle_reactions = []

    def has_heterocycle(smiles):
        """Check if a molecule contains any heterocyclic ring"""
        for ring in HETEROCYCLIC_RINGS:
            if checker.check_ring(ring, smiles):
                if ring not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append(ring)
                return True
        return False

    def is_leaf_node(node):
        """Check if a node is a leaf node (no children)"""
        return len(node.get("children", [])) == 0

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_starting_materials, convergent_heterocycle_reactions, findings_json

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if this is a starting material with a heterocycle
            if node.get("in_stock", False) and has_heterocycle(mol_smiles):
                heterocycle_starting_materials.append(mol_smiles)
            elif is_leaf_node(node) and has_heterocycle(mol_smiles):
                # Leaf nodes without children are also considered starting materials
                heterocycle_starting_materials.append(mol_smiles)

        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    return

                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check if this is a convergent reaction (multiple reactants)
                if len(reactants) >= 2:
                    # Check if at least one reactant contains a heterocycle
                    heterocycle_reactants = []
                    for r in reactants:
                        if has_heterocycle(r):
                            heterocycle_reactants.append(r)

                    # Check if the product contains a heterocycle
                    product_has_heterocycle = has_heterocycle(product_part)

                    # Check if this is a heterocycle-forming reaction
                    is_heterocycle_forming = False
                    for rxn_type in HETEROCYCLE_FORMING_REACTIONS:
                        if checker.check_reaction(rxn_type, rsmi):
                            is_heterocycle_forming = True
                            if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                            break

                    # If we have heterocycle reactants and either the product has a heterocycle
                    # or this is a known heterocycle-forming reaction
                    if len(heterocycle_reactants) > 0 and (
                        product_has_heterocycle or is_heterocycle_forming
                    ):
                        convergent_heterocycle_reactions.append(rsmi)
            except Exception:
                pass

        # Continue traversal
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "mol"
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have multiple heterocycle starting materials and at least one convergent reaction
    result = len(heterocycle_starting_materials) >= 2 and len(convergent_heterocycle_reactions) >= 1

    # Populate structural constraints based on the result
    if len(heterocycle_starting_materials) >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "starting_material_with_heterocycle",
                "operator": ">=",
                "value": 2
            }
        })
    if len(convergent_heterocycle_reactions) >= 1:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "convergent_heterocyclic_reaction",
                "operator": ">=",
                "value": 1
            }
        })

    return result, findings_json

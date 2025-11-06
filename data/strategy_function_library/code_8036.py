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


HETEROCYCLE_RING_NAMES = [
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
    "phenoxazine",
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
    "pyrazole",
    "Paal-Knorr pyrrole",
    "Fischer indole",
    "oxadiazole",
    "benzofuran",
    "benzothiophene",
    "indole",
    "Friedlaender chinoline",
    "triaryl-imidazole",
    "imidazole",
    "Pictet-Spengler",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Identifies routes employing a 'convergent heterocycle assembly' strategy. A route is flagged if it contains BOTH: (1) at least one convergent reaction step that couples two or more fragments already containing rings from the `HETEROCYCLE_RING_NAMES` list, and (2) at least one reaction step that forms a new heterocycle. Heterocycle formation is identified either by the appearance of a new ring scaffold (in multi-component reactions only) from the `HETEROCYCLE_RING_NAMES` list or by matching a named reaction from the `HETEROCYCLE_FORMING_REACTIONS` list.
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

    convergent_steps = 0
    heterocycle_formations = 0

    def dfs_traverse(node, depth=0):
        nonlocal convergent_steps, heterocycle_formations, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            if len(reactants_smiles) > 1:
                heterocycle_reactants = 0
                reactant_heterocycles = set()

                for reactant in reactants_smiles:
                    reactant_has_heterocycle = False
                    for ring_name in HETEROCYCLE_RING_NAMES:
                        if checker.check_ring(ring_name, reactant):
                            reactant_has_heterocycle = True
                            if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                            reactant_heterocycles.add(ring_name)

                    if reactant_has_heterocycle:
                        heterocycle_reactants += 1

                product_has_heterocycle = False
                product_heterocycles = set()

                for ring_name in HETEROCYCLE_RING_NAMES:
                    if checker.check_ring(ring_name, product_smiles):
                        product_has_heterocycle = True
                        if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                        product_heterocycles.add(ring_name)

                if heterocycle_reactants >= 2 and product_has_heterocycle:
                    convergent_steps += 1
                    print(f"Found convergent heterocycle step: {rsmi}")

                new_heterocycles = product_heterocycles - reactant_heterocycles
                if new_heterocycles:
                    heterocycle_formations += 1
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                    print(
                        f"Found heterocycle formation (new type): {rsmi}, new heterocycles: {new_heterocycles}"
                    )

            for rxn_type in HETEROCYCLE_FORMING_REACTIONS:
                if checker.check_reaction(rxn_type, rsmi):
                    heterocycle_formations += 1
                    if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    print(f"Found heterocycle-forming reaction {rxn_type}: {rsmi}")
                    break

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = convergent_steps > 0 and heterocycle_formations > 0
    if result:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "convergent_heterocycle_coupling",
                    "heterocycle_formation"
                ],
                "description": "The route must contain at least one instance of a 'convergent_heterocycle_coupling' event AND at least one instance of a 'heterocycle_formation' event. A 'convergent_heterocycle_coupling' is a multi-component reaction where at least two reactants contain a heterocycle. A 'heterocycle_formation' is a reaction that either matches a named heterocycle-forming reaction or results in a new heterocycle scaffold in the product."
            }
        })

    print(
        f"Convergent heterocycle synthesis detected: {result} (convergent steps: {convergent_steps}, heterocycle formations: {heterocycle_formations})"
    )
    return result, findings_json

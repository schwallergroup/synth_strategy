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

HETEROCYCLE_RINGS = [
    "isoxazole", "pyrazole", "triazole", "tetrazole", "oxazole", "thiazole",
    "imidazole", "pyrimidine", "pyridine", "furan", "thiophene", "pyrrole",
    "oxadiazole", "thiadiazole", "benzoxazole", "benzothiazole", "benzimidazole",
    "indole", "quinoline", "isoquinoline", "purine", "piperidine", "piperazine",
    "morpholine", "thiomorpholine", "pyrrolidine", "azetidine", "aziridine",
    "oxirane", "thiirane", "dioxane", "dioxolane", "dithiolane",
    "thiazolidine", "oxazolidine", "pyrazine", "pyridazine", "carbazole",
    "acridine", "benzotriazole", "indazole", "pteridin", "dibenzofuran",
    "dibenzothiophene",
]

HETEROCYCLE_FORMATION_REACTIONS = [
    "{tetrazole_terminal}", "{tetrazole_connect_regioisomere_1}",
    "{tetrazole_connect_regioisomere_2}", "{Huisgen_Cu-catalyzed_1,4-subst}",
    "{1,2,4-triazole_acetohydrazide}", "{pyrazole}", "{oxadiazole}",
    "{benzimidazole_derivatives_carboxylic-acid/ester}",
    "{benzimidazole_derivatives_aldehyde}", "{benzothiazole}",
    "{benzoxazole_arom-aldehyde}", "{benzoxazole_carboxylic-acid}",
    "{thiazole}", "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition", "Pyrazole formation",
    "{Niementowski_quinazoline}", "{Fischer indole}", "{indole}",
    "{Friedlaender chinoline}", "{benzofuran}", "{benzothiophene}",
    "Formation of NOS Heterocycles", "Paal-Knorr pyrrole synthesis",
    "{Paal-Knorr pyrrole}", "{triaryl-imidazole}", "{imidazole}",
    "Benzothiazole formation from aldehyde",
    "Benzothiazole formation from acyl halide",
    "Benzothiazole formation from ester/carboxylic acid",
    "Benzoxazole formation from aldehyde",
    "Benzoxazole formation from acyl halide",
    "Benzoxazole formation from ester/carboxylic acid",
    "Benzoxazole formation (intramolecular)",
    "Benzimidazole formation from aldehyde",
    "Benzimidazole formation from acyl halide",
    "Benzimidazole formation from ester/carboxylic acid",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Azide-nitrile click cycloaddition to tetrazole",
    "Azide-nitrile click cycloaddition to triazole",
    "Intramolecular amination of azidobiphenyls (heterocycle formation)",
    "Intramolecular amination (heterocycle formation)",
    "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects syntheses that construct at least two heterocyclic rings. This is
    triggered by reactions that form a new heterocycle, identified either by
    matching a known named reaction from HETEROCYCLE_FORMATION_REACTIONS or by
    the appearance of a new ring structure from the HETEROCYCLE_RINGS list.
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

    heterocycle_formations = []

    heterocycle_rings = HETEROCYCLE_RINGS
    heterocycle_formation_rxns = HETEROCYCLE_FORMATION_REACTIONS

    def dfs_traverse(node, depth=0, path=None):
        if path is None:
            path = []

        nonlocal heterocycle_formations, findings_json

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                reactants_heterocycles = {}
                for reactant in reactants:
                    for ring in heterocycle_rings:
                        if checker.check_ring(ring, reactant):
                            if ring not in reactants_heterocycles:
                                reactants_heterocycles[ring] = 0
                            reactants_heterocycles[ring] += 1

                product_heterocycles = {}
                for ring in heterocycle_rings:
                    if checker.check_ring(ring, product):
                        product_heterocycles[ring] = 1

                formed_heterocycles = []
                for ring in product_heterocycles:
                    if (
                        ring not in reactants_heterocycles
                        or product_heterocycles[ring] > reactants_heterocycles[ring]
                    ):
                        formed_heterocycles.append(ring)
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)

                if formed_heterocycles:
                    for rxn_type in heterocycle_formation_rxns:
                        if checker.check_reaction(rxn_type, rsmi):
                            heterocycle_formations.append((depth, formed_heterocycles, rxn_type))
                            if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                            break
                    else:
                        heterocycle_formations.append((depth, formed_heterocycles, "generic"))
                        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth if node["type"] == "reaction" else depth + 1
            dfs_traverse(child, new_depth, path + [node])

    dfs_traverse(route)

    result = len(heterocycle_formations) >= 2

    if result:
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "heterocycle_formation",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json

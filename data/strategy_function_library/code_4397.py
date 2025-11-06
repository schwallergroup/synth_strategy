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


NITROGEN_HETEROCYCLES_OF_INTEREST = [
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
    "indole",
    "quinoline",
    "isoquinoline",
    "benzimidazole",
    "benzoxazole",
    "benzothiazole",
    "piperidine",
    "piperazine",
    "morpholine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the formation of specific nitrogen-containing heterocycles. The function identifies these formations by checking for a list of relevant named reactions (e.g., Paal-Knorr, Fischer indole) or by confirming the presence of a heterocycle from a defined list (e.g., pyrrole, pyridine, indole) in the product that is absent in the reactants.
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

    has_heterocycle_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_heterocycle_formation, findings_json

        if has_heterocycle_formation:
            return

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                if not product_smiles or not reactants_smiles:
                    return

                # Check for specific heterocycle formation reactions
                for rxn_name in [
                    "Paal-Knorr pyrrole synthesis",
                    "benzimidazole_derivatives_carboxylic-acid/ester",
                    "benzimidazole_derivatives_aldehyde",
                    "benzothiazole",
                    "benzoxazole_arom-aldehyde",
                    "benzoxazole_carboxylic-acid",
                    "thiazole",
                    "tetrazole_terminal",
                    "1,2,4-triazole_acetohydrazide",
                    "1,2,4-triazole_carboxylic-acid/ester",
                    "pyrazole",
                    "Fischer indole",
                    "indole",
                    "oxadiazole",
                    "imidazole",
                ]:
                    if checker.check_reaction(rxn_name, rsmi):
                        has_heterocycle_formation = True
                        if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                        return

                # Check for the formation of specific nitrogen-containing heterocycles
                for ring_name in NITROGEN_HETEROCYCLES_OF_INTEREST:
                    if checker.check_ring(ring_name, product_smiles):
                        reactant_has_ring = any(checker.check_ring(ring_name, r) for r in reactants_smiles)
                        if not reactant_has_ring:
                            has_heterocycle_formation = True
                            if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                            return
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return has_heterocycle_formation, findings_json

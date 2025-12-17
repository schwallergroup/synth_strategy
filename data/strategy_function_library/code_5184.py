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


# Refactored module-level constants
HETEROCYCLES_OF_INTEREST = [
    "quinoline", "isoquinoline", "pyridine", "pyrimidine", "pyrazine",
    "pyridazine", "indole", "benzimidazole", "benzoxazole", "benzothiazole",
    "furan", "thiophene", "pyrrole", "imidazole", "oxazole", "thiazole",
    "triazole", "tetrazole", "purine",
]

FG_REACTIONS_OF_INTEREST = [
    "Oxidation of aldehydes to carboxylic acids",
    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
    "Reduction of aldehydes and ketones to alcohols",
    "Reduction of carboxylic acid to primary alcohol",
    "Reduction of ester to primary alcohol",
    "Reduction of nitrile to amine",
    "Reduction of nitro groups to amines",
    "Esterification of Carboxylic Acids",
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Alcohol to chloride_SOCl2",
    "Alcohol to azide",
    "Azide to amine reduction (Staudinger)",
    "Boc amine protection",
    "Boc amine deprotection",
    "Methylation",
    "Alkylation of amines",
    "Williamson Ether Synthesis",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route follows a functional group interconversion cascade
    on a preserved heterocyclic core.
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

    fg_transformations = 0
    heterocycle_core_preserved = True
    active_heterocycle = None
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal fg_transformations, heterocycle_core_preserved, active_heterocycle, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                try:
                    product_has_core = False
                    reactant_has_core = False

                    if active_heterocycle is None:
                        for reactant in reactants:
                            for heterocycle in HETEROCYCLES_OF_INTEREST:
                                if checker.check_ring(heterocycle, reactant):
                                    active_heterocycle = heterocycle
                                    reactant_has_core = True
                                    if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                        findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                                    break
                            if active_heterocycle:
                                break
                    else:
                        for reactant in reactants:
                            if checker.check_ring(active_heterocycle, reactant):
                                reactant_has_core = True
                                if active_heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                                    findings_json["atomic_checks"]["ring_systems"].append(active_heterocycle)
                                break

                    if active_heterocycle:
                        product_has_core = checker.check_ring(active_heterocycle, product)

                    if reactant_has_core and not product_has_core:
                        heterocycle_core_preserved = False
                        # Record ring destruction
                        if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")

                    if reactant_has_core and product_has_core:
                        for reaction_type in FG_REACTIONS_OF_INTEREST:
                            if checker.check_reaction(reaction_type, rsmi):
                                fg_transformations += 1
                                if reaction_type not in findings_json["atomic_checks"]["named_reactions"]:
                                    findings_json["atomic_checks"]["named_reactions"].append(reaction_type)
                                break
                except Exception:
                    pass

        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if fg_transformations >= 3 and heterocycle_core_preserved and active_heterocycle is not None:
        result = True
        # Record structural constraints if the overall condition is met
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "functional_group_interconversion_on_preserved_core",
                "operator": ">=",
                "value": 3
            }
        })
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "ring_destruction",
                "scope": "identified_heterocyclic_core"
            }
        })
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "any_heterocycle_of_interest",
                "operator": ">=",
                "value": 1
            }
        })

    return result, findings_json
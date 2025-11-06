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


def is_piperazine_functionalization(reaction, depth, max_depth) -> Tuple[bool, Dict]:
    """
    Checks if a reaction is a functionalization of an existing piperazine scaffold, meaning the moiety is present in both a reactant and the product.
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
    result = False

    # A functionalization requires the piperazine scaffold to be present before and after the reaction.
    product_has_piperazine = checker.has_group(reaction['product'], 'piperazine')
    if product_has_piperazine:
        findings_json["atomic_checks"]["ring_systems"].append("piperazine")

    reactant_has_piperazine = False
    for reactant in reaction['reactants']:
        if checker.has_group(reactant, 'piperazine'):
            reactant_has_piperazine = True
            findings_json["atomic_checks"]["ring_systems"].append("piperazine")
            break

    result = product_has_piperazine and reactant_has_piperazine

    if result:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "piperazine_in_product",
                    "piperazine_in_reactant"
                ],
                "scope": "single_step"
            }
        })

    return result, findings_json
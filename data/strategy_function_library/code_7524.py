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


# Refactored lists for enumeration
AZIDE_FORMATION_REACTIONS = [
    "Formation of Azides from halogens",
    "Formation of Azides from boronic acids",
    "Alcohol to azide",
    "Amine to azide",
]

AZIDE_TRANSFORMATION_REACTIONS = [
    "Azide to amine reduction (Staudinger)",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Azide-nitrile click cycloaddition to tetrazole",
    "Azide-nitrile click cycloaddition to triazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where an azide is used as a synthetic handle. This is confirmed by identifying at least one reaction that forms an azide and at least one subsequent reaction that transforms it into a different nitrogen-containing group.
    The function specifically checks for known azide formation reactions (e.g., from halides, boronic acids, alcohols, amines) and known azide transformations (e.g., Staudinger reduction, Huisgen cycloadditions, azide-nitrile cycloadditions) using curated lists of reaction names.
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

    # Track molecules containing azides and their transformations
    azide_molecules = set()
    azide_transformations = []

    def dfs_traverse(node, depth=0):
        nonlocal azide_molecules, azide_transformations, findings_json
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            product_has_azide = checker.check_fg("Azide", product)
            if product_has_azide:
                if "Azide" not in findings_json["atomic_checks"]["functional_groups"]:
                    findings_json["atomic_checks"]["functional_groups"].append("Azide")

            # Check for azide formation (retrosynthetically: some precursor -> azide)
            if product_has_azide:
                for rxn in AZIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        azide_molecules.add(product)
                        if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        break

            # Check for azide transformation (retrosynthetically: other N group <- azide)
            # Product must not be an azide, but must contain nitrogen, and a reactant must be an azide.
            if not product_has_azide and "N" in product:
                reactant_has_azide = False
                for reactant in reactants:
                    if checker.check_fg("Azide", reactant):
                        reactant_has_azide = True
                        if "Azide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Azide")
                        break

                if reactant_has_azide:
                    # Check if this is a known azide transformation reaction
                    for rxn in AZIDE_TRANSFORMATION_REACTIONS:
                        if checker.check_reaction(rxn, rsmi):
                            # Find the specific azide reactant
                            for reactant in reactants:
                                if checker.check_fg("Azide", reactant):
                                    azide_transformations.append((reactant, product, rsmi))
                                    if rxn not in findings_json["atomic_checks"]["named_reactions"]:
                                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                                    break
                            break

        # Traverse children (deeper in retrosynthesis = earlier in forward synthesis)
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # We want both azide formation and subsequent transformation
    has_azide_formation = len(azide_molecules) > 0
    has_azide_transformation = len(azide_transformations) > 0

    result = has_azide_formation and has_azide_transformation

    if result:
        # Add the structural constraint if both conditions are met
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "azide_formation_reaction",
                    "azide_transformation_reaction"
                ],
                "description": "The route must contain at least one reaction from the azide formation list (e.g., 'Formation of Azides from halogens') and at least one reaction from the azide transformation list (e.g., 'Azide to amine reduction (Staudinger)')."
            }
        })

    return result, findings_json

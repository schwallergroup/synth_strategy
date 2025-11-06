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


CYCLIZATION_REACTION_NAMES = [
    "Formation of NOS Heterocycles",
    "Paal-Knorr pyrrole synthesis",
    "Intramolecular amination (heterocycle formation)",
    "Intramolecular amination of azidobiphenyls (heterocycle formation)",
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
    "Pictet-Spengler",
    "Diels-Alder",
    "Diels-Alder (ON bond)",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Pyrazole formation",
    "Michael-induced ring closure from hydrazone",
    "Michael-induced ring closure from diazoalkane",
    "[3+2]-cycloaddition of hydrazone and alkyne",
    "[3+2]-cycloaddition of hydrazone and alkene",
    "[3+2]-cycloaddition of diazoalkane and alkyne",
    "[3+2]-cycloaddition of diazoalkane and alkene",
    "[3+2]-cycloaddition of diazoalkane and alpha-alkyne",
    "[3+2]-cycloaddition of diazoalkane and alpha-alkene",
    "Retro-Diels-Alder from oxazole",
    "Pauson-Khand reaction",
    "Azide-nitrile click cycloaddition to tetrazole",
    "Azide-nitrile click cycloaddition to triazole",
    "Alkyne-imine cycloaddition",
    "A3 coupling to imidazoles",
    "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
    "Intramolecular transesterification/Lactone formation",
    "{Pictet-Spengler}",
    "{benzimidazole_derivatives_carboxylic-acid/ester}",
    "{benzimidazole_derivatives_aldehyde}",
    "{benzothiazole}",
    "{benzoxazole_arom-aldehyde}",
    "{benzoxazole_carboxylic-acid}",
    "{thiazole}",
    "{Niementowski_quinazoline}",
    "{tetrazole_terminal}",
    "{tetrazole_connect_regioisomere_1}",
    "{tetrazole_connect_regioisomere_2}",
    "{Huisgen_Cu-catalyzed_1,4-subst}",
    "{Huisgen_Ru-catalyzed_1,5_subst}",
    "{Huisgen_disubst-alkyne}",
    "{1,2,4-triazole_acetohydrazide}",
    "{1,2,4-triazole_carboxylic-acid/ester}",
    "{3-nitrile-pyridine}",
    "{spiro-chromanone}",
    "{pyrazole}",
    "{phthalazinone}",
    "{Paal-Knorr pyrrole}",
    "{triaryl-imidazole}",
    "{Fischer indole}",
    "{Friedlaender chinoline}",
    "{benzofuran}",
    "{benzothiophene}",
    "{indole}",
    "{oxadiazole}",
    "{imidazole}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a late-stage (final step) cyclization reaction. This is identified by first checking if the reaction matches a comprehensive list of known cyclization reaction names. If no match is found, it performs a rigorous analysis using atom-mapping numbers to confirm that at least one atom moves from an acyclic environment in the reactants to a cyclic one in the product.
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

    has_late_stage_cyclization = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_cyclization, findings_json

        if node["type"] == "reaction" and depth <= 1:
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Analyzing late-stage reaction at depth {depth}: {rsmi}")

                # Check 1: Use a predefined list of known cyclization reaction types
                for rxn in CYCLIZATION_REACTION_NAMES:
                    if checker.check_reaction(rxn, rsmi):
                        print(f"Detected known cyclization reaction type: {rxn}")
                        has_late_stage_cyclization = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        # If a known cyclization reaction is found, we can consider the 'ring_formation' constraint met.
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "ring_formation",
                                "position": "last_stage"
                            }
                        })
                        return

                # Check 2: Analyze atom mappings to rigorously detect ring formation
                if "[" in rsmi and ":" in rsmi:  # Check if atom mapping exists
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    if not all(reactant_mols) or not product_mol:
                        print("Warning: Could not parse all molecules")
                        return

                    product_ring_info = product_mol.GetRingInfo()
                    product_atom_rings = {}

                    for atom in product_mol.GetAtoms():
                        atom_idx = atom.GetIdx()
                        if product_ring_info.NumAtomRings(atom_idx) > 0:
                            if atom.HasProp("molAtomMapNumber"):
                                map_num = atom.GetProp("molAtomMapNumber")
                                product_atom_rings[map_num] = atom_idx

                    for reactant_mol in reactant_mols:
                        if not reactant_mol:
                            continue
                        reactant_ring_info = reactant_mol.GetRingInfo()
                        for atom in reactant_mol.GetAtoms():
                            if atom.HasProp("molAtomMapNumber"):
                                map_num = atom.GetProp("molAtomMapNumber")
                                atom_idx = atom.GetIdx()
                                if (
                                    map_num in product_atom_rings
                                    and reactant_ring_info.NumAtomRings(atom_idx) == 0
                                ):
                                    print(
                                        f"Atom with mapping {map_num} forms a new ring in the product"
                                    )
                                    has_late_stage_cyclization = True
                                    # Record the 'ring_formation' atomic check
                                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                                    # Record the structural constraint
                                    findings_json["structural_constraints"].append({
                                        "type": "positional",
                                        "details": {
                                            "target": "ring_formation",
                                            "position": "last_stage"
                                        }
                                    })
                                    return
            except Exception as e:
                print(f"Error analyzing reaction for cyclization: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    print(f"Late-stage cyclization strategy detected: {has_late_stage_cyclization}")
    return has_late_stage_cyclization, findings_json

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


FG_TYPES_OF_INTEREST = [
    "Primary alcohol",
    "Secondary alcohol",
    "Tertiary alcohol",
    "Aromatic alcohol",
    "Aldehyde",
    "Ketone",
    "Carboxylic acid",
    "Ester",
    "Acyl halide",
    "Anhydride",
    "Nitrile",
    "Primary amine",
    "Secondary amine",
    "Tertiary amine",
    "Primary amide",
    "Secondary amide",
    "Tertiary amide",
    "Nitro group",
    "Primary halide",
    "Secondary halide",
    "Tertiary halide",
    "Aromatic halide",
    "Boronic acid",
    "Boronic ester",
    "Alkyne",
    "Azide",
    "Isocyanate",
    "Sulfonamide",
    "Sulfonic acid",
    "Sulfonate",
    "Thiol",
    "Phosphate ester",
]

FG_INTERCONVERSION_REACTIONS = [
    "Oxidation of aldehydes to carboxylic acids",
    "Oxidation of alcohols to carboxylic acids",
    "Oxidation of alkene to carboxylic acid",
    "Oxidation of ketone to carboxylic acid",
    "Oxidation of nitrile to carboxylic acid",
    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
    "Oxidation of primary alcohols",
    "Reduction of aldehydes and ketones to alcohols",
    "Reduction of carboxylic acid to primary alcohol",
    "Reduction of ester to primary alcohol",
    "Reduction of nitrile to amine",
    "Reduction of primary amides to amines",
    "Reduction of secondary amides to amines",
    "Reduction of tertiary amides to amines",
    "Esterification of Carboxylic Acids",
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
    "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a cycle of functional group interconversions on the same atom within a synthesis route. This is achieved by first identifying reactions that match a predefined list of common interconversions (`FG_INTERCONVERSION_REACTIONS`). For each matched reaction, it then determines the specific functional groups involved (from the `FG_TYPES_OF_INTEREST` list) on mapped atoms. Finally, it analyzes the sequence of transformations for each atom across the entire route to detect futile cycles (e.g., A -> B -> A).
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

    # Track functional group transformations by atom mapping
    atom_transformations = {}

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # This means it's a 'chemical' node
            next_depth = depth + 1

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            # Parse reactants and products
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_part.split(".")]
            product_mol = Chem.MolFromSmiles(products_part)

            if None in reactant_mols or product_mol is None:
                print(f"Warning: Could not parse molecules in reaction: {rsmi}")
                for child in node.get("children", []):
                    dfs_traverse(child, next_depth)
                return

            # Check for functional groups in reactants and products
            reactant_fgs = {}
            product_fgs = {}

            # Check functional groups in reactants
            for r_mol_idx, r_mol in enumerate(reactant_mols):
                r_smiles = Chem.MolToSmiles(r_mol)
                for fg in FG_TYPES_OF_INTEREST:
                    if checker.check_fg(fg, r_smiles):
                        if fg not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append(fg)
                        fg_atoms = checker.get_fg_atom_indices(fg, r_smiles)
                        for atom_indices in fg_atoms:
                            for atom in atom_indices:
                                for a in r_mol.GetAtoms():
                                    if a.GetIdx() == atom:
                                        if a.HasProp("molAtomMapNumber"):
                                            map_num = int(a.GetProp("molAtomMapNumber"))
                                            if map_num not in reactant_fgs:
                                                reactant_fgs[map_num] = []
                                            reactant_fgs[map_num].append(fg)
                                        break

            # Check functional groups in product
            p_smiles = Chem.MolToSmiles(product_mol)
            for fg in FG_TYPES_OF_INTEREST:
                if checker.check_fg(fg, p_smiles):
                    if fg not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append(fg)
                    fg_atoms = checker.get_fg_atom_indices(fg, p_smiles)
                    for atom_indices in fg_atoms:
                        for atom in atom_indices:
                            for a in product_mol.GetAtoms():
                                if a.GetIdx() == atom:
                                    if a.HasProp("molAtomMapNumber"):
                                        map_num = int(a.GetProp("molAtomMapNumber"))
                                        if map_num not in product_fgs:
                                            product_fgs[map_num] = []
                                        product_fgs[map_num].append(fg)
                                    break

            # Check if this reaction is a functional group interconversion
            is_fg_interconversion = False
            for rxn_type in FG_INTERCONVERSION_REACTIONS:
                if checker.check_reaction(rxn_type, rsmi):
                    if rxn_type not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_type)
                    is_fg_interconversion = True
                    print(f"Detected reaction type: {rxn_type}")
                    break

            if is_fg_interconversion:
                for map_num in set(reactant_fgs.keys()).union(set(product_fgs.keys())):
                    r_fgs = reactant_fgs.get(map_num, [])
                    p_fgs = product_fgs.get(map_num, [])

                    if r_fgs and p_fgs and set(r_fgs) != set(p_fgs):
                        for r_fg in r_fgs:
                            for p_fg in p_fgs:
                                if r_fg != p_fg:
                                    transformation = f"{r_fg}_to_{p_fg}"
                                    print(
                                        f"Detected transformation: {transformation} on atom {map_num} at depth {depth}"
                                    )

                                    if map_num not in atom_transformations:
                                        atom_transformations[map_num] = []
                                    atom_transformations[map_num].append((transformation, depth))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check for cycles in transformations
    has_cycle = False
    for atom_map, transformations in atom_transformations.items():
        if len(transformations) >= 2:
            # Record the first structural constraint if met
            constraint_1 = {
                "type": "count",
                "details": {
                    "description": "At least two functional group interconversion reactions must occur on the same atom to enable a cycle.",
                    "target": "functional_group_interconversion_on_same_atom",
                    "operator": ">=",
                    "value": 2
                }
            }
            if constraint_1 not in findings_json["structural_constraints"]:
                findings_json["structural_constraints"].append(constraint_1)

            sorted_trans = sorted(transformations, key=lambda x: x[1], reverse=True)
            print(f"Atom {atom_map} transformations: {sorted_trans}")

            fg_sequence = []
            for trans in sorted_trans:
                trans_type = trans[0]
                start_fg = trans_type.split("_to_")[0]
                end_fg = trans_type.split("_to_")[1]

                if not fg_sequence:
                    fg_sequence.append(start_fg)
                elif fg_sequence[-1] != start_fg:
                    fg_sequence.append(start_fg)

                fg_sequence.append(end_fg)

            print(f"Functional group sequence for atom {atom_map}: {fg_sequence}")

            for i, fg in enumerate(fg_sequence):
                if fg in fg_sequence[:i]:
                    print(
                        f"Found cycle for atom {atom_map}: functional group {fg} appears multiple times"
                    )
                    print(f"Full sequence: {' -> '.join(fg_sequence)}")
                    has_cycle = True
                    # Record the second structural constraint if met
                    constraint_2 = {
                        "type": "sequence",
                        "details": {
                            "description": "The ordered sequence of functional groups on a single atom, generated by interconversions, must contain a repeated functional group, indicating a futile cycle.",
                            "target": "sequence_of_functional_groups_on_a_single_atom",
                            "property": "contains_duplicates"
                        }
                    }
                    if constraint_2 not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(constraint_2)
                    break

            if not has_cycle:
                for i in range(len(sorted_trans) - 1):
                    trans1 = sorted_trans[i][0]
                    trans2 = sorted_trans[i + 1][0]

                    start_fg1 = trans1.split("_to_")[0]
                    end_fg1 = trans1.split("_to_")[1]

                    start_fg2 = trans2.split("_to_")[0]
                    end_fg2 = trans2.split("_to_")[1]

                    if end_fg1 == start_fg2 and end_fg2 == start_fg1:
                        print(
                            f"Found direct cycle for atom {atom_map}: {start_fg1} → {end_fg1} → {start_fg1}"
                        )
                        has_cycle = True
                        # Record the second structural constraint if met
                        constraint_2 = {
                            "type": "sequence",
                            "details": {
                                "description": "The ordered sequence of functional groups on a single atom, generated by interconversions, must contain a repeated functional group, indicating a futile cycle.",
                                "target": "sequence_of_functional_groups_on_a_single_atom",
                                "property": "contains_duplicates"
                            }
                        }
                        if constraint_2 not in findings_json["structural_constraints"]:
                            findings_json["structural_constraints"].append(constraint_2)
                        break

    print(f"Has functional group cycle: {has_cycle}")
    return has_cycle, findings_json

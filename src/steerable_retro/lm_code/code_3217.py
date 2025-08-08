#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
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

root_data = "/home/dparm/steerable_retro/data"

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


def main(route):
    """
    This function detects if the synthesis uses a late-stage ring formation strategy,
    particularly focusing on the final steps (depth <= 2).
    """
    # Track if we found late-stage ring formation
    has_late_stage_ring_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_ring_formation

        # Check reactions in the late stage (final step and up to 2 steps before)
        if node["type"] == "reaction" and depth <= 2:
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")
                product = product_part

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                product_mol = Chem.MolFromSmiles(product) if product else None
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

                if not product_mol or not reactant_mols:
                    print("Could not parse product or reactants")
                    return

                # Count rings in product and reactants
                product_ring_count = product_mol.GetRingInfo().NumRings()
                reactant_ring_counts = [r.GetRingInfo().NumRings() for r in reactant_mols if r]
                max_reactant_ring_count = max(reactant_ring_counts, default=0)

                print(f"Product ring count: {product_ring_count}")
                print(f"Reactant ring counts: {reactant_ring_counts}")
                print(f"Max reactant ring count: {max_reactant_ring_count}")

                # Check if the reaction is a known ring-forming reaction type
                ring_forming_reaction = False
                ring_forming_reaction_types = [
                    "Formation of NOS Heterocycles",
                    "Paal-Knorr pyrrole synthesis",
                    "Diels-Alder",
                    "Diels-Alder (ON bond)",
                    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
                    "Huisgen 1,3 dipolar cycloaddition",
                    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
                    "Pyrazole formation",
                    "A3 coupling to imidazoles",
                    "Alkyne-imine cycloaddition",
                    "Azide-nitrile click cycloaddition to tetrazole",
                    "Azide-nitrile click cycloaddition to triazole",
                    "Michael-induced ring closure from hydrazone",
                    "Michael-induced ring closure from diazoalkane",
                    "[3+2]-cycloaddition of hydrazone and alkyne",
                    "[3+2]-cycloaddition of hydrazone and alkene",
                    "[3+2]-cycloaddition of diazoalkane and alkyne",
                    "[3+2]-cycloaddition of diazoalkane and alkene",
                    "[3+2]-cycloaddition of diazoalkane and alpha-alkyne",
                    "[3+2]-cycloaddition of diazoalkane and alpha-alkene",
                    "Intramolecular amination of azidobiphenyls (heterocycle formation)",
                    "Intramolecular amination (heterocycle formation)",
                    "Pauson-Khand reaction",
                    "Pictet-Spengler",
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
                    "Huisgen_disubst-alkyne",
                    "1,2,4-triazole_acetohydrazide",
                    "1,2,4-triazole_carboxylic-acid/ester",
                    "3-nitrile-pyridine",
                    "spiro-chromanone",
                    "pyrazole",
                    "phthalazinone",
                    "Paal-Knorr pyrrole",
                    "triaryl-imidazole",
                    "Fischer indole",
                    "Friedlaender chinoline",
                    "benzofuran",
                    "benzothiophene",
                    "indole",
                    "oxadiazole",
                    "imidazole",
                    "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
                ]

                for rxn_type in ring_forming_reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected ring-forming reaction: {rxn_type}")
                        ring_forming_reaction = True
                        break

                # Check for new ring types in product that weren't in reactants
                common_ring_types = [
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
                    "benzothiophene",
                    "oxathiolane",
                    "dioxathiolane",
                    "thiazolidine",
                    "oxazolidine",
                    "isoxazole",
                    "isothiazole",
                    "oxadiazole",
                    "thiadiazole",
                    "cyclopropane",
                    "cyclobutane",
                    "cyclopentane",
                    "cyclohexane",
                    "cycloheptane",
                    "cyclooctane",
                    "benzene",
                    "naphthalene",
                    "anthracene",
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

                product_rings = set()
                for ring_type in common_ring_types:
                    if checker.check_ring(ring_type, product):
                        product_rings.add(ring_type)
                        print(f"Found {ring_type} in product")

                reactant_rings = set()
                for reactant in reactants:
                    for ring_type in common_ring_types:
                        if checker.check_ring(ring_type, reactant):
                            reactant_rings.add(ring_type)
                            print(f"Found {ring_type} in reactant")

                new_rings = product_rings - reactant_rings
                if new_rings:
                    print(f"New ring types formed: {new_rings}")

                # Check for intramolecular ring formation
                intramolecular_ring_formation = False
                if len(reactants) == 1 and product_ring_count > max_reactant_ring_count:
                    print("Detected intramolecular ring formation")
                    intramolecular_ring_formation = True

                # Check for ring fusion
                # This is a simplistic approach - if product has rings and more atoms in rings than reactants
                product_ring_atoms = set()
                for ring in product_mol.GetRingInfo().AtomRings():
                    product_ring_atoms.update(ring)

                reactant_ring_atoms_count = 0
                for r_mol in reactant_mols:
                    for ring in r_mol.GetRingInfo().AtomRings():
                        reactant_ring_atoms_count += len(ring)

                ring_fusion = (
                    len(product_ring_atoms) > reactant_ring_atoms_count
                    and product_ring_count <= max_reactant_ring_count
                )
                if ring_fusion:
                    print("Detected potential ring fusion")

                # Determine if this is a late-stage ring formation
                if (
                    (product_ring_count > max_reactant_ring_count)
                    or ring_forming_reaction
                    or new_rings
                    or intramolecular_ring_formation
                    or ring_fusion
                ):
                    print(f"Found late-stage ring formation at depth {depth}")
                    has_late_stage_ring_formation = True
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    if route["type"] == "mol":
        dfs_traverse(route)
    else:
        print("Root is not a molecule node")

    return has_late_stage_ring_formation

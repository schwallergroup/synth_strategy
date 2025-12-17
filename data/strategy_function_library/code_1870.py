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


HETEROCYCLES_OF_INTEREST = [
    "quinoline", "isoquinoline", "pyridine", "pyrimidine", "pyrazine",
    "pyridazine", "triazole", "tetrazole", "benzimidazole", "benzoxazole",
    "benzothiazole", "indole", "purine", "quinazoline", "pteridin",
    "imidazole", "oxazole", "thiazole", "furan", "thiophene", "pyrrole",
]

LINKER_FGS = [
    "Primary halide", "Secondary halide", "Tertiary halide", "Ether",
    "Primary alcohol", "Secondary alcohol", "Tertiary alcohol", "Primary amine",
    "Secondary amine", "Tertiary amine", "Nitrile", "Carboxylic acid",
    "Ester", "Amide", "Phenol", "Aromatic halide",
]

LINKER_INSTALLATION_REACTIONS = [
    "Williamson Ether Synthesis", "Williamson Ether Synthesis (intra to epoxy)",
    "N-alkylation of primary amines with alkyl halides",
    "N-alkylation of secondary amines with alkyl halides", "S-alkylation of thiols",
    "S-alkylation of thiols (ethyl)", "S-alkylation of thiols with alcohols",
    "S-alkylation of thiols with alcohols (ethyl)", "Esterification of Carboxylic Acids",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", "Suzuki coupling with boronic acids",
    "Suzuki coupling with boronic esters", "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
]

HETEROCYCLE_FORMATION_REACTIONS = [
    "Formation of NOS Heterocycles", "Paal-Knorr pyrrole synthesis",
    "benzimidazole_derivatives_carboxylic-acid/ester", "benzimidazole_derivatives_aldehyde",
    "benzothiazole", "benzoxazole_arom-aldehyde", "benzoxazole_carboxylic-acid",
    "thiazole", "Niementowski_quinazoline", "tetrazole_terminal",
    "tetrazole_connect_regioisomere_1", "tetrazole_connect_regioisomere_2",
    "Huisgen_Cu-catalyzed_1,4-subst", "Huisgen_Ru-catalyzed_1,5_subst",
    "1,2,4-triazole_acetohydrazide", "1,2,4-triazole_carboxylic-acid/ester",
    "pyrazole", "Paal-Knorr pyrrole", "Fischer indole", "Friedlaender chinoline",
    "benzofuran", "benzothiophene", "indole", "oxadiazole",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition", "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition", "Pyrazole formation",
]

LINKER_MODIFICATION_REACTIONS = [
    "Oxidation of aldehydes to carboxylic acids", "Oxidation of alcohol to carboxylic acid",
    "Reduction of ester to primary alcohol", "Reduction of ketone to secondary alcohol",
    "Reduction of carboxylic acid to primary alcohol", "Alcohol to azide",
    "Nitrile to amide", "Alcohol to ether", "Ester saponification (methyl deprotection)",
    "Ester saponification (alkyl deprotection)",
    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
    "Reduction of primary amides to amines", "Reduction of secondary amides to amines",
    "Reduction of tertiary amides to amines", "Reduction of nitrile to amine",
    "Primary amine to fluoride", "Primary amine to chloride", "Primary amine to bromide",
    "Primary amine to iodide", "Alcohol to chloride_sulfonyl chloride",
    "Alcohol to chloride_SOCl2", "Alcohol to triflate conversion", "Appel reaction",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a multi-step strategy involving: 1. Installation of a linker functional group via a defined set of reactions. 2. Formation of a specific heterocycle on the linker-containing intermediate. 3. Optional, subsequent modification of the linker via a defined set of reactions. The strategy is only flagged if the linker installation occurs synthetically before the heterocycle formation. The specific linkers, heterocycles, and reactions are defined in module-level lists.
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

    # Track if we found the key features
    found_linker_installation = False
    found_heterocycle_formation = False
    found_linker_modification = False

    # Track the depth at which each feature was found
    linker_installation_depth = -1
    heterocycle_formation_depth = -1
    linker_modification_depth = -1

    # Track the molecules with installed linkers and heterocycles
    molecules_with_linkers = {}  # molecule SMILES -> set of linker FGs
    molecules_with_heterocycles = {}  # molecule SMILES -> set of heterocycles

    # Track the synthesis path
    synthesis_path = []

    def dfs_traverse(node, depth=0, path=None):
        nonlocal found_linker_installation, found_heterocycle_formation, found_linker_modification
        nonlocal linker_installation_depth, heterocycle_formation_depth, linker_modification_depth
        nonlocal findings_json

        if path is None:
            path = []

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                is_linker_installation = False
                for rxn in LINKER_INSTALLATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        is_linker_installation = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        break

                if is_linker_installation:
                    found_linker_installation = True
                    linker_installation_depth = depth
                    linker_fgs_in_product = set()
                    for fg in LINKER_FGS:
                        if checker.check_fg(fg, product_smiles):
                            linker_fgs_in_product.add(fg)
                            findings_json["atomic_checks"]["functional_groups"].append(fg)
                    if linker_fgs_in_product:
                        molecules_with_linkers[product_smiles] = linker_fgs_in_product

                is_heterocycle_formation = False
                for rxn in HETEROCYCLE_FORMATION_REACTIONS:
                    if checker.check_reaction(rxn, rsmi):
                        is_heterocycle_formation = True
                        findings_json["atomic_checks"]["named_reactions"].append(rxn)
                        break

                if not is_heterocycle_formation:
                    new_heterocycles = set()
                    for ring in HETEROCYCLES_OF_INTEREST:
                        if checker.check_ring(ring, product_smiles) and not any(checker.check_ring(ring, r) for r in reactants_smiles):
                            new_heterocycles.add(ring)
                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                    is_heterocycle_formation = len(new_heterocycles) > 0

                if is_heterocycle_formation:
                    reactants_with_linkers = [r for r in reactants_smiles if r in molecules_with_linkers]

                    if reactants_with_linkers:
                        product_has_linker = any(checker.check_fg(fg, product_smiles) for fg in LINKER_FGS)

                        if product_has_linker:
                            found_heterocycle_formation = True
                            heterocycle_formation_depth = depth
                            heterocycles_in_product = set()
                            for ring in HETEROCYCLES_OF_INTEREST:
                                if checker.check_ring(ring, product_smiles):
                                    heterocycles_in_product.add(ring)
                                    if ring not in findings_json["atomic_checks"]["ring_systems"]:
                                        findings_json["atomic_checks"]["ring_systems"].append(ring)
                            if heterocycles_in_product:
                                molecules_with_heterocycles[product_smiles] = heterocycles_in_product

                if found_heterocycle_formation:
                    reactants_with_heterocycles = [r for r in reactants_smiles if r in molecules_with_heterocycles]
                    if not reactants_with_heterocycles:
                        for r in reactants_smiles:
                            if any(checker.check_ring(ring, r) for ring in HETEROCYCLES_OF_INTEREST):
                                reactants_with_heterocycles.append(r)
                                break

                    if reactants_with_heterocycles:
                        is_linker_modification = False
                        for rxn in LINKER_MODIFICATION_REACTIONS:
                            if checker.check_reaction(rxn, rsmi):
                                is_linker_modification = True
                                findings_json["atomic_checks"]["named_reactions"].append(rxn)
                                break

                        product_has_heterocycle = any(checker.check_ring(ring, product_smiles) for ring in HETEROCYCLES_OF_INTEREST)

                        if is_linker_modification and product_has_heterocycle:
                            found_linker_modification = True
                            linker_modification_depth = depth
                            synthesis_path.append(
                                {"type": "linker_modification", "depth": depth, "rsmi": rsmi}
                            )
            except Exception:
                pass

        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # This means the current node is a chemical node
                new_depth = depth + 1
            dfs_traverse(child, new_depth, path + [node])

    dfs_traverse(route)

    result = False
    if found_linker_installation and found_heterocycle_formation:
        # Add co-occurrence constraint
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "linker_installation",
                    "heterocycle_formation"
                ],
                "description": "The strategy requires at least one linker installation reaction and one heterocycle formation reaction to occur in the synthesis route."
            }
        })

        if linker_installation_depth > heterocycle_formation_depth:
            result = True
            # Add sequence constraint
            findings_json["structural_constraints"].append({
                "type": "sequence",
                "details": {
                    "before": "linker_installation",
                    "after": "heterocycle_formation",
                    "description": "The linker installation step must occur synthetically before the heterocycle formation step."
                }
            })

    # Deduplicate findings_json lists
    for key in findings_json["atomic_checks"]:
        findings_json["atomic_checks"][key] = list(set(findings_json["atomic_checks"][key]))

    return result, findings_json

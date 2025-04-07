#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    This function detects a synthetic strategy involving multiple late-stage
    functional group modifications (depths 0-2).
    """
    late_stage_modifications = 0

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_modifications

        if node["type"] == "reaction" and depth <= 3:  # Extended to depth 3
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for various functional group modifications using reaction types

                # Methylation reactions
                if (
                    checker.check_reaction("N-methylation", rsmi)
                    or checker.check_reaction("O-methylation", rsmi)
                    or checker.check_reaction("S-methylation", rsmi)
                    or checker.check_reaction("C-methylation", rsmi)
                    or checker.check_reaction("Methylation", rsmi)
                    or checker.check_reaction("Methylation with MeI_primary", rsmi)
                    or checker.check_reaction("Methylation with MeI_secondary", rsmi)
                    or checker.check_reaction("Methylation with MeI_tertiary", rsmi)
                    or checker.check_reaction("Methylation with MeI_aryl", rsmi)
                    or checker.check_reaction("Methylation with MeI_SH", rsmi)
                    or checker.check_reaction("Methylation with DMS", rsmi)
                    or checker.check_reaction("Methylation of OH with DMS", rsmi)
                    or checker.check_reaction("Methylation with DMC", rsmi)
                    or checker.check_reaction("DMS COOH methylation", rsmi)
                    or checker.check_reaction("DMS Amine methylation", rsmi)
                    or checker.check_reaction("Eschweiler-Clarke Primary Amine Methylation", rsmi)
                    or checker.check_reaction("Eschweiler-Clarke Secondary Amine Methylation", rsmi)
                    or checker.check_reaction(
                        "Reductive methylation of primary amine with formaldehyde", rsmi
                    )
                    or checker.check_reaction("Parnes methylation", rsmi)
                ):
                    late_stage_modifications += 1
                    print(f"Late-stage methylation detected at depth {depth}")

                # Halogenation/dehalogenation reactions
                elif (
                    checker.check_reaction("Aromatic fluorination", rsmi)
                    or checker.check_reaction("Aromatic chlorination", rsmi)
                    or checker.check_reaction("Aromatic bromination", rsmi)
                    or checker.check_reaction("Aromatic iodination", rsmi)
                    or checker.check_reaction("Chlorination", rsmi)
                    or checker.check_reaction("Fluorination", rsmi)
                    or checker.check_reaction("Iodination", rsmi)
                    or checker.check_reaction("Bromination", rsmi)
                    or checker.check_reaction("Dehalogenation", rsmi)
                    or checker.check_reaction("Aromatic dehalogenation", rsmi)
                    or checker.check_reaction("Finkelstein reaction", rsmi)
                    or checker.check_reaction("Primary amine to fluoride", rsmi)
                    or checker.check_reaction("Primary amine to chloride", rsmi)
                    or checker.check_reaction("Primary amine to bromide", rsmi)
                    or checker.check_reaction("Primary amine to iodide", rsmi)
                ):
                    late_stage_modifications += 1
                    print(f"Late-stage halogenation/dehalogenation detected at depth {depth}")

                # Oxidation/reduction reactions
                elif (
                    checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi)
                    or checker.check_reaction(
                        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
                    )
                    or checker.check_reaction("Oxidation of ketone to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of alcohol to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of nitrile to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of amide to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of alkene to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of boronic acids", rsmi)
                    or checker.check_reaction("Oxidation of boronic esters", rsmi)
                    or checker.check_reaction("Reduction of ester to primary alcohol", rsmi)
                    or checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi)
                    or checker.check_reaction(
                        "Reduction of carboxylic acid to primary alcohol", rsmi
                    )
                    or checker.check_reaction(
                        "Reduction of aldehydes and ketones to alcohols", rsmi
                    )
                    or checker.check_reaction("Reduction of nitro groups to amines", rsmi)
                    or checker.check_reaction("Reduction of primary amides to amines", rsmi)
                    or checker.check_reaction("Reduction of secondary amides to amines", rsmi)
                    or checker.check_reaction("Reduction of tertiary amides to amines", rsmi)
                    or checker.check_reaction("Reduction of nitrile to amine", rsmi)
                ):
                    late_stage_modifications += 1
                    print(f"Late-stage oxidation/reduction detected at depth {depth}")

                # Protection/deprotection reactions
                elif (
                    checker.check_reaction("Alcohol protection with silyl ethers", rsmi)
                    or checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                    or checker.check_reaction(
                        "Alcohol deprotection from silyl ethers (double)", rsmi
                    )
                    or checker.check_reaction("Alcohol deprotection from silyl ethers (diol)", rsmi)
                    or checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                    or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                    or checker.check_reaction("Boc amine protection", rsmi)
                    or checker.check_reaction("Boc amine protection explicit", rsmi)
                    or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                    or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                    or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                    or checker.check_reaction("Boc amine protection of primary amine", rsmi)
                    or checker.check_reaction("Deprotection of carboxylic acid", rsmi)
                    or checker.check_reaction("Protection of carboxylic acid", rsmi)
                    or checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                    or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                    or checker.check_reaction("Hydroxyl benzyl deprotection", rsmi)
                    or checker.check_reaction("Carboxyl benzyl deprotection", rsmi)
                    or checker.check_reaction("COOH ethyl deprotection", rsmi)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                    or checker.check_reaction("TMS deprotection from alkyne", rsmi)
                    or checker.check_reaction("N-glutarimide deprotection", rsmi)
                    or checker.check_reaction("Phthalimide deprotection", rsmi)
                ):
                    late_stage_modifications += 1
                    print(f"Late-stage protection/deprotection detected at depth {depth}")

                # Esterification/amidation reactions
                elif (
                    checker.check_reaction("Schotten-Baumann to ester", rsmi)
                    or checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    or checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                    )
                    or checker.check_reaction(
                        "O-alkylation of carboxylic acids with diazo compounds", rsmi
                    )
                    or checker.check_reaction("O-alkylation of amides with diazo compounds", rsmi)
                    or checker.check_reaction("Aminolysis of esters", rsmi)
                    or checker.check_reaction("Formation of Sulfonic Esters", rsmi)
                    or checker.check_reaction(
                        "Formation of Sulfonic Esters on TMS protected alcohol", rsmi
                    )
                    or checker.check_reaction("Oxidative esterification of primary alcohols", rsmi)
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction("Acylation of secondary amines with anhydrides", rsmi)
                    or checker.check_reaction("Acylation of secondary amines", rsmi)
                    or checker.check_reaction("Acylation of primary amines", rsmi)
                    or checker.check_reaction("Acyl chloride with ammonia to amide", rsmi)
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                    or checker.check_reaction("Acyl chloride with primary amine to imide", rsmi)
                    or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                    or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with ammonia to amide", rsmi)
                    or checker.check_reaction("Ester with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                    or checker.check_reaction("Nitrile to amide", rsmi)
                    or checker.check_reaction("Carboxylic acid to amide conversion", rsmi)
                    or checker.check_reaction("Acetic anhydride and alcohol to ester", rsmi)
                ):
                    late_stage_modifications += 1
                    print(f"Late-stage esterification/amidation detected at depth {depth}")

                # Alkylation/arylation reactions
                elif (
                    checker.check_reaction("Williamson Ether Synthesis", rsmi)
                    or checker.check_reaction("Williamson Ether Synthesis (intra to epoxy)", rsmi)
                    or checker.check_reaction("S-alkylation of thiols", rsmi)
                    or checker.check_reaction("S-alkylation of thiols (ethyl)", rsmi)
                    or checker.check_reaction("S-alkylation of thiols with alcohols", rsmi)
                    or checker.check_reaction("S-alkylation of thiols with alcohols (ethyl)", rsmi)
                    or checker.check_reaction("Alkylation of amines", rsmi)
                    or checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction("Friedel-Crafts alkylation", rsmi)
                    or checker.check_reaction("Friedel-Crafts alkylation with halide", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                    or checker.check_reaction("Suzuki coupling with sulfonic esters", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                    or checker.check_reaction("Negishi coupling", rsmi)
                    or checker.check_reaction(
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                    )
                    or checker.check_reaction("Stille reaction_vinyl", rsmi)
                    or checker.check_reaction("Stille reaction_aryl", rsmi)
                    or checker.check_reaction("Stille reaction_benzyl", rsmi)
                    or checker.check_reaction("Stille reaction_allyl", rsmi)
                    or checker.check_reaction("Stille reaction_vinyl OTf", rsmi)
                    or checker.check_reaction("Stille reaction_aryl OTf", rsmi)
                    or checker.check_reaction("Stille reaction_benzyl OTf", rsmi)
                    or checker.check_reaction("Stille reaction_allyl OTf", rsmi)
                    or checker.check_reaction("Stille reaction_other", rsmi)
                    or checker.check_reaction("Stille reaction_other OTf", rsmi)
                    or checker.check_reaction("Hiyama-Denmark Coupling", rsmi)
                    or checker.check_reaction("Kumada cross-coupling", rsmi)
                    or checker.check_reaction("Aryllithium cross-coupling", rsmi)
                ):
                    late_stage_modifications += 1
                    print(f"Late-stage alkylation/arylation detected at depth {depth}")

                # Other functional group transformations
                elif (
                    checker.check_reaction("Formation of Azides from halogens", rsmi)
                    or checker.check_reaction("Formation of Azides from boronic acids", rsmi)
                    or checker.check_reaction("Azide to amine reduction (Staudinger)", rsmi)
                    or checker.check_reaction("Alcohol to azide", rsmi)
                    or checker.check_reaction("Amine to azide", rsmi)
                    or checker.check_reaction("Alcohol to chloride_sulfonyl chloride", rsmi)
                    or checker.check_reaction("Alcohol to chloride_SOCl2", rsmi)
                    or checker.check_reaction("Alcohol to chloride_CHCl3", rsmi)
                    or checker.check_reaction("Alcohol to chloride_CH2Cl2", rsmi)
                    or checker.check_reaction("Alcohol to chloride_PCl5_ortho", rsmi)
                    or checker.check_reaction("Alcohol to chloride_POCl3_ortho", rsmi)
                    or checker.check_reaction("Alcohol to chloride_POCl3_para", rsmi)
                    or checker.check_reaction("Alcohol to chloride_POCl3", rsmi)
                    or checker.check_reaction("Alcohol to chloride_HCl", rsmi)
                    or checker.check_reaction("Alcohol to chloride_Salt", rsmi)
                    or checker.check_reaction("Alcohol to chloride_Other", rsmi)
                    or checker.check_reaction("Alcohol to triflate conversion", rsmi)
                    or checker.check_reaction("Alcohol to ether", rsmi)
                    or checker.check_reaction("Cleavage of methoxy ethers to alcohols", rsmi)
                    or checker.check_reaction("Cleavage of alkoxy ethers to alcohols", rsmi)
                    or checker.check_reaction("Ether cleavage to primary alcohol", rsmi)
                    or checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                    )
                    or checker.check_reaction("Sulfanyl to sulfinyl_peroxide", rsmi)
                    or checker.check_reaction("Sulfanyl to sulfinyl_H2O2", rsmi)
                    or checker.check_reaction("Sulfanyl to sulfinyl_H2O", rsmi)
                    or checker.check_reaction("Sulfanyl to sulfinyl_SO3-", rsmi)
                    or checker.check_reaction("Sulfanyl to sulfinyl_sulfonyl", rsmi)
                    or checker.check_reaction("Sulfanyl to sulfinyl_MeOH", rsmi)
                    or checker.check_reaction("Sulfanyl to sulfinyl_COO", rsmi)
                    or checker.check_reaction("Sulfanyl to sulfinyl", rsmi)
                    or checker.check_reaction("PBr3 and alcohol to alkyl bromide", rsmi)
                    or checker.check_reaction("Nitrile and hydrogen peroxide to amide", rsmi)
                    or checker.check_reaction("Appel reaction", rsmi)
                    or checker.check_reaction("Carboxylic acid to carboxylate", rsmi)
                    or checker.check_reaction("Ester to carboxylate", rsmi)
                    or checker.check_reaction("Carboxylate to carboxylic acid", rsmi)
                    or checker.check_reaction(
                        "Aldehyde and ketone to alpha,beta-unsaturated carbonyl", rsmi
                    )
                    or checker.check_reaction("Amine and thiophosgene to isothiocyanate", rsmi)
                    or checker.check_reaction("Hydrazine synthesis from amine", rsmi)
                    or checker.check_reaction("Aromatic sulfonyl chlorination", rsmi)
                    or checker.check_reaction("Ring opening of epoxide with amine", rsmi)
                ):
                    late_stage_modifications += 1
                    print(f"Late-stage functional group transformation detected at depth {depth}")

                # If no specific reaction type was detected, check for functional group changes
                else:
                    try:
                        # Check for functional group changes by comparing reactants and product
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                        product_mol = Chem.MolFromSmiles(product_smiles)

                        # List of functional groups to check
                        fg_list = [
                            "Primary alcohol",
                            "Secondary alcohol",
                            "Tertiary alcohol",
                            "Primary amine",
                            "Secondary amine",
                            "Tertiary amine",
                            "Carboxylic acid",
                            "Ester",
                            "Amide",
                            "Nitrile",
                            "Aldehyde",
                            "Ketone",
                            "Alkyne",
                            "Alkene",
                            "Aromatic halide",
                            "Primary halide",
                            "Secondary halide",
                            "Tertiary halide",
                            "Phenol",
                            "Ether",
                            "Nitro group",
                            "Azide",
                            "Sulfonamide",
                            "Sulfone",
                            "Thiol",
                            "Sulfoxide",
                            "Enamine",
                            "Hydrazone",
                            "Oxime",
                            "Acetal/Ketal",
                            "Triflate",
                            "Tosylate",
                            "Mesylate",
                            "Anhydride",
                            "Urea",
                            "Thiourea",
                            "Isocyanate",
                            "Isothiocyanate",
                        ]

                        # Check if any functional group is present in reactants but not in product or vice versa
                        fg_changed = False
                        for fg in fg_list:
                            # Check if FG is in any reactant
                            fg_in_reactants = any(
                                checker.check_fg(fg, Chem.MolToSmiles(r))
                                for r in reactant_mols
                                if r is not None
                            )
                            # Check if FG is in product
                            fg_in_product = (
                                checker.check_fg(fg, Chem.MolToSmiles(product_mol))
                                if product_mol is not None
                                else False
                            )

                            if fg_in_reactants != fg_in_product:
                                fg_changed = True
                                print(f"Late-stage {fg} modification detected at depth {depth}")
                                late_stage_modifications += 1
                                break
                    except Exception as e:
                        print(f"Error checking functional groups at depth {depth}: {e}")

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if at least 1 late-stage modification is found (changed from 2)
    strategy_present = late_stage_modifications >= 1

    if strategy_present:
        print(
            f"Late-stage functional group modification strategy detected with {late_stage_modifications} modifications"
        )
    else:
        print("Late-stage functional group modification strategy not detected")

    return strategy_present

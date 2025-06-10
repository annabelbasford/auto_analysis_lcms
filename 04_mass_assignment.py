import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
from pathlib import Path
import rdkit 
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem import Draw
from pyopenms import EmpiricalFormula, FineIsotopePatternGenerator
from pyopenms import *
import csv
from math import isclose
#from tqdm import tqdm
import json
from pathlib import Path
from typing import List
from matplotlib.offsetbox import OffsetImage, AnchoredOffsetbox, HPacker

'''
import faulthandler
faulthandler.enable()

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
'''
def _get_precursor_formula(smiles: str) -> EmpiricalFormula:
    precursor = CalcMolFormula(rdkit.Chem.MolFromSmiles(smiles))
    for i in precursor:
        if i == '-':
            precursor = precursor.replace('-', '') # to avoid error of -ve charge in EmpricialFormula
    precursor_formula = EmpiricalFormula(precursor)
    precursor_formula.setCharge(0) # to avoid error of -ve charge in EmpricialFormula
    return precursor_formula

def poc_find_solutions(max_precursors,topicity_a, topicity_b):
    solutions = []
    for a in range(1, max_precursors +1):  # set number of topicimer A
        for b in range(0, max_precursors +1):  # set number topicimer B1
            if b * topicity_b < a * topicity_a:
                max_imines = b * topicity_b
            else:
                max_imines = a * topicity_a
            for x in range(a+b-1,max_imines+1): #set connectivity from linear connectivity up to maximum number of imines
                if a + b < max_precursors:
                    line = a,b,x #create tuple of values
                    solutions.append(line)
    solutions.sort(key=lambda x: x[0])  # Sort solutions by aldehydes
    return solutions

def poc_find_solutions_ternary(max_precursors,topicity_a, topicity_b1, topicity_b2):
    solutions = []
    for a in range(1, max_precursors +1):  # set number of topicimer A
        for b1 in range(0, max_precursors +1):  # set number topicimer B1
            for b2 in range(0, max_precursors +1): # set number of topicimer B2
                if b1 * topicity_b1 + b2 * topicity_b2 < a * topicity_a:
                    max_imines = b1 * topicity_b1 + b2 * topicity_b2
                else:
                    max_imines = a * topicity_a
                for x in range(a+b1+b2-1,max_imines+1): #set connectivity from linear connectivity up to maximum number of imines
                    if a + b1 + b2 < max_precursors:
                        line = a,b1,b2,x #create tuple of values
                        solutions.append(line)
    solutions.sort(key=lambda x: x[0])  # Sort solutions by aldehydes
    return solutions

def poc_calc_formulas(solutions, prec_a_smiles, prec_b_smiles,max_charge):
    prec_a_formula= _get_precursor_formula(prec_a_smiles)
    prec_b1_formula = _get_precursor_formula(prec_b_smiles)
    formulas = []
    names = []
    H2O = EmpiricalFormula('H2O')
    for solution in solutions:
        number_of_prec_a = solution[0]
        number_of_prec_b1 = solution[1]
        number_of_imines = solution[2]
        compound_formula= EmpiricalFormula(prec_a_formula) #start generating Empirical formula with openpyms with one precursor a
        for i in range(0,number_of_prec_a-1): #for each component start adding on fragments
                compound_formula += prec_a_formula
        for i in range(0,number_of_prec_b1):
                compound_formula += prec_b1_formula
        for i in range(0,number_of_imines):
                compound_formula -= H2O
        for charge in range(1,max_charge+1):
            compound_formula += EmpiricalFormula('H')
            gen_name = 'X'+str(number_of_prec_a)+'_Y'+str(number_of_prec_b1)+'_Bonds'+str(number_of_imines) + '_Charge'+str(charge) #generate name for dictionary key
            names.append(gen_name)
            comp_properties = {'formula':compound_formula,'charge': charge}
            formulas.append(comp_properties)
    data = dict(zip(names,formulas))
    return data

def poc_calc_formulas_ternary(solutions, prec_a_smiles, prec_b1_smiles, prec_b2_smiles,max_charge):
    prec_a_formula= _get_precursor_formula(prec_a_smiles)
    prec_b1_formula = _get_precursor_formula(prec_b1_smiles)
    prec_b2_formula = _get_precursor_formula(prec_b2_smiles) 
    formulas = []
    names = []
    H2O = EmpiricalFormula('H2O')
    for solution in solutions:
        number_of_prec_a = solution[0]
        number_of_prec_b1 = solution[1]
        number_of_prec_b2 = solution[2]
        number_of_imines = solution[3]
        compound_formula=prec_a_formula #start generating Empirical formula with openpyms with one precursor a
        for i in range(0,number_of_prec_a-1): #for each component start adding on fragments
                compound_formula = compound_formula + prec_a_formula
        for i in range(0,number_of_prec_b1):
                compound_formula = compound_formula + prec_b1_formula
        for i in range(0,number_of_prec_b2):
                compound_formula = compound_formula + prec_b2_formula
        for i in range(0,number_of_imines):
                compound_formula = compound_formula - H2O
        for charge in range(1,max_charge+1):
            compound_formula = compound_formula + EmpiricalFormula('H')
            gen_name = 'X'+str(number_of_prec_a)+'_Y'+str(number_of_prec_b1)+'_Z'+str(number_of_prec_b2)+'_Bonds'+str(number_of_imines) + '_Charge'+str(charge) #generate name for dictionary key
            names.append(gen_name)
            comp_properties = {'formula':compound_formula,'charge': charge}
            formulas.append(comp_properties) 
    data = dict(zip(names,formulas))
    return data

def poc_calc_formulas_full_cages(solutions, prec_a_smiles, topicity_a, prec_b_smiles, topicity_b, max_charge):
    prec_a_formula= _get_precursor_formula(prec_a_smiles)
    prec_b1_formula = _get_precursor_formula(prec_b_smiles)
    formulas = []
    names = []
    H2O = EmpiricalFormula('H2O')
    for solution in solutions:
        number_of_prec_a = solution[0]
        number_of_prec_b1 = solution[1]
        number_of_imines = solution[2]
        if number_of_prec_a * topicity_a == number_of_prec_b1 * topicity_b == number_of_imines:
            compound_formula= EmpiricalFormula(prec_a_formula) #start generating Empirical formula with openpyms with one precursor a
            for i in range(0,number_of_prec_a-1): #for each component start adding on fragments
                    compound_formula += prec_a_formula
            for i in range(0,number_of_prec_b1):
                    compound_formula += prec_b1_formula
            for i in range(0,number_of_imines):
                    compound_formula -= H2O
            for charge in range(1,max_charge+1):
                compound_formula += EmpiricalFormula('H')
                gen_name = 'X'+str(number_of_prec_a)+'_Y'+str(number_of_prec_b1)+'_Bonds'+str(number_of_imines) + '_Charge'+str(charge) #generate name for dictionary key
                names.append(gen_name)
                comp_properties = {'formula':compound_formula,'charge': charge}
                formulas.append(comp_properties)
    data = dict(zip(names,formulas))
    return data

def get_top_10_isotopes(formula_dict, error=1e-2): # error previously 1e-3
    top_10_isotopes_dict = {}

#    with tqdm(total=len(formula_dict), desc="Calulating isotope peaks") as pbar:
    for key, value in formula_dict.items():
        formula = value['formula']
        charge = value['charge']

        # Generate isotope distribution with the threshold
        isotopes = formula.getIsotopeDistribution(FineIsotopePatternGenerator(error, False))

        # Normalize m/z values by charge
        normalized_isotopes = [(iso.getMZ() / charge, iso.getIntensity()) for iso in isotopes.getContainer()]

        # Extract and sort isotopes by intensity
        normalized_isotopes.sort(key=lambda x: x[1], reverse=True)

        # Get top 10 highest abundance isotopes
        top_10_isotopes = normalized_isotopes[:10]

        # Store the top 10 isotopes in the dictionary
        top_10_isotopes_dict[key] = {
            'predicted_charge': charge,
            'isotopes': top_10_isotopes
        }
#            pbar.update(1)

    return top_10_isotopes_dict

def read_csv(file_path):
    df = pd.read_csv(file_path)

    # Normalize column names to lowercase
    df.columns = df.columns.str.lower()

    # Rename 'm/z' to 'mz' if present
    if 'm/z' in df.columns:
        df = df.rename(columns={'m/z': 'mz'})

    # Now check for expected column combinations
    if 'mz' in df.columns and 'height' in df.columns:
        df = df.rename(columns={'height': 'intensity'})
        mz_col = 'mz'
        intensity_col = 'intensity'
    elif 'mz' in df.columns and 'intensity' in df.columns:
        mz_col = 'mz'
        intensity_col = 'intensity'
    else:
        raise ValueError(f"File {file_path} does not contain required 'mz' or 'm/z' and 'height' or 'intensity' columns.")

    return list(zip(df[mz_col], df[intensity_col]))
'''
def read_csv(file_path):
    df = pd.read_csv(file_path)

    # Normalize column names to lowercase and strip spaces
    df.columns = df.columns.str.strip().str.lower()

    # Show column names for debugging
    print(f"Column names after normalization: {df.columns.tolist()}")

    # Rename 'm/z' to 'mz' explicitly if present (after lowercasing)
    column_renames = {}
    for col in df.columns:
        if col.replace(" ", "") in ["m/z", "mz"]:
            column_renames[col] = "mz"
        elif col == "height":
            column_renames[col] = "intensity"

    df = df.rename(columns=column_renames)

    # Final check
    print(f"Final column names after renaming: {df.columns.tolist()}")
    print(df)
    if 'mz' in df.columns and 'intensity' in df.columns:
        return list(zip(df['mz'], df['intensity']))
    else:
        raise ValueError(f"File {file_path} does not contain required 'mz' or 'm/z' and 'height' or 'intensity' columns.")
'''

def find_matching_isotopes_df(mz_intensity_data, top_10_isotopes_dict, mass_error_ppm=20):
    matching_results = []
    
    # Determine the maximum intensity in the data
    max_intensity = max(intensity for mz, intensity in mz_intensity_data)
    intensity_threshold = 0.01* max_intensity  # 0.5% of the maximum intensity previous 
    
    #with tqdm(total=len(mz_intensity_data), desc="Processing peaks") as pbar:
    for mz, intensity in mz_intensity_data:
        if intensity >= intensity_threshold:  # Filtering peaks with intensity above 5% of max intensity
            for formula, isotopes_info in top_10_isotopes_dict.items():
                predicted_charge = isotopes_info['predicted_charge']
                isotopes = isotopes_info['isotopes']
                for iso_mz, iso_intensity in isotopes:
                    theoretical_mass = iso_mz
                    mass_error = theoretical_mass * mass_error_ppm * 1e-6
                    if isclose(mz, theoretical_mass, abs_tol=mass_error):
                        matching_results.append({
                            'Formula': formula,
                            'found_mz': mz,
                            'theoretical_mz': iso_mz,
                            'found_intensity': intensity,
                            'theoretical_intensity': iso_intensity,
                            'predicted_charge': predicted_charge
                        })
        #pbar.update(1)
    
    # Convert list of dictionaries to DataFrame
    isotopes_df = pd.DataFrame(matching_results)
    return isotopes_df


def calculate_mz_differences(isotopes_df):
    # Sort by 'Formula' and 'found_mz' to ensure correct difference calculation
    df = isotopes_df.sort_values(by=['Formula', 'found_mz'],ascending=False)

    # Create new columns to store the mz differences and charges
    df['mz_difference'] = None
    df['charge'] = None
    df['splitting_error'] = None

    tolerance = 0.05
    
    # Iterate over the DataFrame grouped by 'Formula'
    for formula, group in df.groupby('Formula'):
        previous_mz = None
        for index, row in group.iterrows():
            if previous_mz is not None:
                mz_diff = row['found_mz'] - previous_mz
                df.at[index, 'mz_difference'] = mz_diff
                mz_diff_charge_working = row['predicted_charge']*mz_diff
                if round(abs(mz_diff_charge_working)) == 1: # product of charge and splitting must be equal to 1
                    if abs(1-abs(mz_diff_charge_working)) < tolerance: 
                        df.at[index, 'charge'] = row['predicted_charge']
                        df.at[index, 'splitting_error'] = abs(1-abs(mz_diff_charge_working))*100 # percentage error away from expected mz diff 
            previous_mz = row['found_mz']

    return df

def get_csv_file_names(directory: Path) -> List[Path]:
    return [file for file in directory.iterdir() if file.suffix == '.csv']


def save_json(data, filename):
    """Save data to a JSON file."""
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4)

'''
def main(csv_dir: str, Aldehyde_SMILES, Amine1_SMILES, Amine2_SMILES=None):
    max_precursors = 30

    # Get the data paths and load the data from the reactions json file
    csv_directory = Path(csv_dir)
    csv_file_path_list = get_csv_file_names(csv_directory)

    amine_pattern = Chem.MolFromSmarts('[NH2]')
    aldehyde_pattern = Chem.MolFromSmarts('[CX3H1](=O)[#6]')

    prec_a_smiles = Aldehyde_SMILES
    topicity_a = len(Chem.MolFromSmiles(prec_a_smiles).GetSubstructMatches(aldehyde_pattern))
    prec_b1_smiles = Amine1_SMILES
    topicity_b1 = len(Chem.MolFromSmiles(prec_b1_smiles).GetSubstructMatches(amine_pattern))
    
    if Amine2_SMILES is not None:
        prec_b2_smiles = Amine2_SMILES
        topicity_b2 = len(Chem.MolFromSmiles(prec_b2_smiles).GetSubstructMatches(amine_pattern))

    # get the stems of the csv filenames to match with the codes in reactions_data
    csv_file_stems = {f.stem for f in csv_file_path_list}
           
    for n,csv_path in enumerate(csv_file_path_list):
        
        print(f"Processing {csv_path.name}")

        if Amine2_SMILES is not None:
            precursor_combinations = poc_find_solutions_ternary(max_precursors, topicity_a, topicity_b1, topicity_b2)
            poc_formula_dict = poc_calc_formulas_ternary(precursor_combinations, prec_a_smiles, prec_b1_smiles, prec_b2_smiles, 4)
        else:
            precursor_combinations = poc_find_solutions(max_precursors, topicity_a, topicity_b1)
            poc_formula_dict = poc_calc_formulas(precursor_combinations, prec_a_smiles, prec_b_smiles, 4)
        isotopes_dict = get_top_10_isotopes(poc_formula_dict)
        mz_intensity_data = read_csv(csv_path)
        matching_peaks_df = find_matching_isotopes_df(mz_intensity_data, isotopes_dict)
        if matching_peaks_df.empty:
            print(f"Skipping {csv_path.name}: No matching peaks found.")
            continue
        result_df = calculate_mz_differences(matching_peaks_df)
        result_filtered_df = result_df[result_df['predicted_charge'] == result_df['charge']]

        output_file = csv_path.with_name(f"{csv_path.stem}_output.csv")
        result_filtered_df.to_csv(output_file, index=False)
        print(f"Saved: {output_file}")

        # Reload, sort and group data
        df = pd.read_csv(csv_path)
        if 'm/z' in df.columns:
            df = df.rename(columns={'m/z': 'mz'})
        df2 = result_filtered_df.sort_values(by=["found_intensity"], ascending=False)
        highest_intensity_peaks = df2.groupby("Formula")["found_intensity"].idxmax()

        # Create plot
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.bar(df["mz"], df["Intensity"], width=0.01, edgecolor="red", color="red", alpha=0.3)

        # Add peak labels
        for idx, row in df2.loc[highest_intensity_peaks].iterrows():
            ax.bar(row["found_mz"], row["found_intensity"], width=0.01, edgecolor="black", color="black")
            ax.text(row["found_mz"], row["found_intensity"], row["Formula"], 
                    ha="center", va="bottom", fontsize=8, clip_on=True)

        # Save processed output
        df2.to_csv(str(Path(csv_path).with_name(Path(csv_path).stem+'_processed').with_suffix('.csv')), index=False)

        # get images of the precursor structures using rdkit
        smiles_list = [prec_a_smiles, prec_b1_smiles, prec_b2_smiles]
        mol_images = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                img = Draw.MolToImage(mol, size=(300, 300))
                mol_images.append(OffsetImage(img, zoom=0.6))

        # Add insets to top-right corner
        if mol_images:
            hbox = HPacker(children=mol_images, align="center", pad=0, sep=10)
            anchored_box = AnchoredOffsetbox(
                loc='upper right', child=hbox, pad=0.5, frameon=False,
                bbox_to_anchor=(1, 1), bbox_transform=ax.transAxes, borderpad=0.3
            )
            ax.add_artist(anchored_box)

        # Final plot formatting
        ax.set_xlim(199.5, 3200.5)
        ax.set_xlabel("m/z")
        ax.set_ylabel("Intensity")
        ax.set_title("Matching Peaks (Highest Intensity)")
        plt.title(str(Path(csv_path).stem))
        plt.tight_layout()
        plt.savefig(str(Path(csv_path).with_suffix('.png')),dpi=300)
        #plt.show()
        plt.close()
'''

def main(csv_dir: str, Aldehyde_SMILES, Amine1_SMILES, Amine2_SMILES=None):
    max_precursors = 12

    # Set up paths
    csv_directory = Path(csv_dir)
    output_dir = csv_directory / "assigned_peaks_data"
    output_dir.mkdir(exist_ok=True)  # Create if it doesn't exist

    csv_file_path_list = get_csv_file_names(csv_directory)

    amine_pattern = Chem.MolFromSmarts('[NH2]')
    aldehyde_pattern = Chem.MolFromSmarts('[CX3H1](=O)[#6]')

    prec_a_smiles = Aldehyde_SMILES
    topicity_a = len(Chem.MolFromSmiles(prec_a_smiles).GetSubstructMatches(aldehyde_pattern))
    prec_b1_smiles = Amine1_SMILES
    topicity_b1 = len(Chem.MolFromSmiles(prec_b1_smiles).GetSubstructMatches(amine_pattern))

    if Amine2_SMILES is not None:
        prec_b2_smiles = Amine2_SMILES
        topicity_b2 = len(Chem.MolFromSmiles(prec_b2_smiles).GetSubstructMatches(amine_pattern))

    for csv_path in csv_file_path_list:
        print(f"Processing {csv_path.name}")

        if Amine2_SMILES is not None:
            precursor_combinations = poc_find_solutions_ternary(max_precursors, topicity_a, topicity_b1, topicity_b2)
            poc_formula_dict = poc_calc_formulas_ternary(precursor_combinations, prec_a_smiles, prec_b1_smiles, prec_b2_smiles, 4)
        else:
            precursor_combinations = poc_find_solutions(max_precursors, topicity_a, topicity_b1)
            poc_formula_dict = poc_calc_formulas(precursor_combinations, prec_a_smiles, prec_b1_smiles, 4)

        isotopes_dict = get_top_10_isotopes(poc_formula_dict)
        mz_intensity_data = read_csv(csv_path)
        matching_peaks_df = find_matching_isotopes_df(mz_intensity_data, isotopes_dict)

        if matching_peaks_df.empty:
            print(f"Skipping {csv_path.name}: No matching peaks found.")
            continue

        result_df = calculate_mz_differences(matching_peaks_df)
        result_filtered_df = result_df[result_df['predicted_charge'] == result_df['charge']]

        # Save filtered results
        output_file = output_dir / f"{csv_path.stem}_output.csv"
        result_filtered_df.to_csv(output_file, index=False)
        print(f"Saved: {output_file}")

        # Reload and prepare plot data
        df = pd.read_csv(csv_path)
        if 'm/z' in df.columns:
            df = df.rename(columns={'m/z': 'mz'})
        df2 = result_filtered_df.sort_values(by=["found_intensity"], ascending=False)
        highest_intensity_peaks = df2.groupby("Formula")["found_intensity"].idxmax()

        # Plot setup
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.bar(df["mz"], df["Intensity"], width=0.01, edgecolor="red", color="red", alpha=0.3)

        for idx, row in df2.loc[highest_intensity_peaks].iterrows():
            ax.bar(row["found_mz"], row["found_intensity"], width=0.01, edgecolor="black", color="black")
            ax.text(row["found_mz"], row["found_intensity"], row["Formula"],
                    ha="center", va="bottom", fontsize=8, clip_on=True)

        # Add RDKit molecule images
        smiles_list = [prec_a_smiles, prec_b1_smiles] + ([prec_b2_smiles] if Amine2_SMILES else [])
        mol_images = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                img = Draw.MolToImage(mol, size=(300, 300))
                mol_images.append(OffsetImage(img, zoom=0.4))

        if mol_images:
            hbox = HPacker(children=mol_images, align="center", pad=0, sep=10)
            anchored_box = AnchoredOffsetbox(
                loc='upper right', child=hbox, pad=0.5, frameon=False,
                bbox_to_anchor=(1, 1), bbox_transform=ax.transAxes, borderpad=0.3
            )
            ax.add_artist(anchored_box)

        # Format and save plot
        ax.set_xlim(199.5, 3200.5)
        ax.set_xlabel("m/z")
        ax.set_ylabel("Intensity")
        ax.set_title(f"Matching Peaks: {csv_path.stem}")
        plt.tight_layout()
        plt.savefig(output_dir / f"{csv_path.stem}.png", dpi=300)
        plt.close()

        # Save processed data
        df2_file = output_dir / f"{csv_path.stem}_processed.csv"
        df2.to_csv(df2_file, index=False)

# Example usage:
# main("/path/to/ms_csv_dir", "/path/to/reactions_json_dir")
if __name__ == "__main__":
    main('/Users/user/Documents/GitHub/auto_analysis_lcms/peak_mz_refined','O=CC1=CC(C=O)=CC(C=O)=C1','NC1C(N)CCCC1','CC(C)(N)CN')

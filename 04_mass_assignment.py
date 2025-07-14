import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import rdkit 
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem import Draw
from pyopenms import EmpiricalFormula, FineIsotopePatternGenerator
from pyopenms import *
from math import isclose
import json
from pathlib import Path
from typing import List
from matplotlib.offsetbox import OffsetImage, AnchoredOffsetbox, HPacker
import os

def _get_precursor_formula(smiles: str) -> EmpiricalFormula:
    precursor = CalcMolFormula(rdkit.Chem.MolFromSmiles(smiles))
    for i in precursor:
        if i == '-':
            precursor = precursor.replace('-', '') # to avoid error of -ve charge in EmpricialFormula
    precursor_formula = EmpiricalFormula(precursor)
    precursor_formula.setCharge(0) # to avoid error of -ve charge in EmpricialFormula
    return precursor_formula

def poc_find_solutions(max_precursors, topicity_aldehyde, topicity_amine):
    """
    Generate all valid precursor combinations of aldehydes and amines
    based on their topicities (number of reactive sites) and a maximum total number of precursors.

    Returns:
        List of tuples: (num_aldehydes, num_amines, num_imine_bonds)
    """
    solutions = []

    for num_aldehydes in range(1, max_precursors + 1):
        for num_amines in range(0, max_precursors + 1):

            # Determine max number of imine bonds that can form based on topicities
            max_possible_bonds = min(
                num_aldehydes * topicity_aldehyde,
                num_amines * topicity_amine
            )

            # Allow bond counts starting from minimum required for a connected network
            min_bonds = num_aldehydes + num_amines - 1

            for num_bonds in range(min_bonds, max_possible_bonds + 1):
                if num_aldehydes + num_amines < max_precursors:
                    solutions.append((num_aldehydes, num_amines, num_bonds))

    # Optional: sort by number of aldehydes
    solutions.sort(key=lambda combo: combo[0])

    return solutions

def poc_find_solutions_ternary(max_precursors, topicity_aldehyde, topicity_amine1, topicity_amine2):
    """
    Generate all valid precursor combinations involving:
    - one aldehyde type
    - two distinct amine types (amine1, amine2)

    Each precursor has a specified number of reactive sites (topicity).
    The total number of precursors must be less than max_precursors.

    Returns:
        List of tuples: (num_aldehydes, num_amines_1, num_amines_2, num_imine_bonds)
    """
    solutions = []

    for num_aldehydes in range(1, max_precursors + 1):
        for num_amines_1 in range(0, max_precursors + 1):
            for num_amines_2 in range(0, max_precursors + 1):

                # Calculate the maximum number of bonds based on topicities
                total_amine_sites = (num_amines_1 * topicity_amine1) + (num_amines_2 * topicity_amine2)
                total_aldehyde_sites = num_aldehydes * topicity_aldehyde
                max_possible_bonds = min(total_aldehyde_sites, total_amine_sites)

                min_bonds = num_aldehydes + num_amines_1 + num_amines_2 - 1

                for num_bonds in range(min_bonds, max_possible_bonds + 1):
                    if num_aldehydes + num_amines_1 + num_amines_2 < max_precursors:
                        solutions.append((num_aldehydes, num_amines_1, num_amines_2, num_bonds))

    # Optional: sort by number of aldehydes
    solutions.sort(key=lambda combo: combo[0])

    return solutions

def poc_calc_formulas(solutions, smiles_aldehyde, smiles_amine, max_charge):
    """
    Given precursor combinations and SMILES strings, generate formulas
    for all charged species using OpenMS's EmpiricalFormula.

    Returns a dictionary mapping a string key to a dict with 'formula' and 'charge'.
    """
    formula_aldehyde = _get_precursor_formula(smiles_aldehyde)
    formula_amine = _get_precursor_formula(smiles_amine)
    formulas = []
    names = []
    H2O = EmpiricalFormula('H2O') # Represents water loss per imine bond

    for solution in solutions:
        num_aldehydes = solution[0]
        num_amines = solution[1]
        num_bonds = solution[2]
        compound_formula= EmpiricalFormula(formula_aldehyde) #start generating Empirical formula with openpyms with one precursor a
        for i in range(0,num_aldehydes-1): #for each component start adding on fragments
                compound_formula = compound_formula + formula_aldehyde
        for i in range(0,num_amines):
                compound_formula = compound_formula + formula_amine
        for i in range(0,num_bonds):
                compound_formula = compound_formula - H2O # Remove water for each imine
        
        # Add protons for each possible charge state
        for charge in range(1,max_charge+1):
            compound_formula = compound_formula + EmpiricalFormula('H')
            gen_name = 'X'+str(num_aldehydes)+'_Y'+str(num_amines)+'_Bonds'+str(num_bonds) + '_Charge'+str(charge) #generate name for dictionary key
            names.append(gen_name)
            comp_properties = {'formula':compound_formula,'charge': charge}
            formulas.append(comp_properties)
    data = dict(zip(names,formulas))
    return data

def poc_calc_formulas_ternary(solutions, smiles_aldehyde, smiles_amine1, smiles_amine2, max_charge):
    """
    Given precursor combinations and SMILES strings, generate formulas
    for all charged species using OpenMS's EmpiricalFormula.

    Returns a dictionary mapping a string key to a dict with 'formula' and 'charge'.
    """
    formula_aldehyde = _get_precursor_formula(smiles_aldehyde)
    formula_amine1 = _get_precursor_formula(smiles_amine1)
    formula_amine2 = _get_precursor_formula(smiles_amine2) 
    formulas = []
    names = []
    H2O = EmpiricalFormula('H2O')
    for solution in solutions:
        num_aldehydes = solution[0]
        num_amines_1 = solution[1]
        num_amines_2 = solution[2]
        num_bonds = solution[3]
        compound_formula=formula_aldehyde #start generating Empirical formula with openpyms with one precursor a
        for i in range(0,num_aldehydes-1): #for each component start adding on fragments
                compound_formula = compound_formula + formula_aldehyde
        for i in range(0,num_amines_1):
                compound_formula = compound_formula + formula_amine1
        for i in range(0,num_amines_2):
                compound_formula = compound_formula + formula_amine2
        for i in range(0,num_bonds):
                compound_formula = compound_formula - H2O
        for charge in range(1,max_charge+1):
            compound_formula = compound_formula + EmpiricalFormula('H')
            gen_name = 'X'+str(num_aldehydes)+'_Y'+str(num_amines_1)+'_Z'+str(num_amines_2)+'_Bonds'+str(num_bonds) + '_Charge'+str(charge) #generate name for dictionary key
            names.append(gen_name)
            comp_properties = {'formula':compound_formula,'charge': charge}
            formulas.append(comp_properties) 
    data = dict(zip(names,formulas))
    return data

def get_top_10_isotopes(formula_dict, error=1e-2): # error previously 1e-3
    top_10_isotopes_dict = {}

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
    """
    Reads a CSV file and standardizes column names for m/z and intensity.
    Accepts both 'm/z' and 'height' as valid intensity representations.
    """
    df = pd.read_csv(file_path)
    df.columns = df.columns.str.lower()

    # Rename 'm/z' to 'mz' if present
    if 'm/z' in df.columns:
        df = df.rename(columns={'m/z': 'mz'})

    if 'Intensity' in df.columns:
        df = df.rename(columns={'Intensity': 'intensity'})

    # Determine final column names for processing
    if 'mz' in df.columns and 'height' in df.columns:
        df = df.rename(columns={'height': 'intensity'})
    elif not ('mz' in df.columns and 'intensity' in df.columns):
        raise ValueError(f"File {file_path} does not contain required columns.")

    return list(zip(df['mz'], df['intensity']))

def find_matching_isotopes_df(mz_intensity_data, top_10_isotopes_dict, mass_error_ppm=20):
    """
    Match observed m/z-intensity peaks to predicted isotopic m/z values within a ppm tolerance.

    Parameters:
        mz_intensity_data: List of tuples (mz, intensity)
        top_10_isotopes_dict: Dictionary of predicted isotopes
        mass_error_ppm: Mass accuracy threshold in parts per million

    Returns:
        A pandas DataFrame of matching peaks and their metadata.
    """
    matching_results = []
    
    # Determine the maximum intensity in the data
    max_intensity = max(intensity for mz, intensity in mz_intensity_data)
    intensity_threshold = 0.01* max_intensity  # 1% of the maximum intensity previous 
    
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
    """
    For each formula group, calculate the m/z difference between consecutive peaks,
    infer charge from spacing, and compute the error from ideal spacing (1 m/z unit per charge).

    Returns:
        DataFrame with additional columns: mz_difference, charge, and splitting_error.
    """
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
    """
    Retrieve all CSV files in the given directory.

    Returns:
        List of file paths with .csv extension.
    """
    return [file for file in directory.iterdir() if file.suffix == '.csv']

def load_experiment_data(json_path: Path) -> List[dict]:
    with open(json_path, 'r') as f:
        return json.load(f)

def main():
    max_precursors = 12 # Maximum number of precursors to consider in the reaction, increasing this will increase the number of possible combinations and thus the time it takes to run the script

    # Set up paths
    parent_folder = '/home/abasford/projects/auto_analysis_lcms/'
    input_folder = '/home/abasford/projects/auto_analysis_lcms/peaks_mz_refined/' # Folder containing the refined CSV files for each peak for each reaction sample
    output_folder = os.path.join(parent_folder, 'assigned_peaks_data') # Folder to save the output files with assigned peaks and plots
    os.makedirs(output_folder, exist_ok=True)

    # read in JSON from the input_folder which maps the raections, number of components with their SMILES
    reactions = load_experiment_data(os.path.join(parent_folder, 'reaction_planner.json'))

    # Get all CSV files in the input folder that end with '_refined.csv'
    csv_file_path_list = [f for f in os.listdir(input_folder) if f.endswith('_refined.csv')]
    print(f"Found {len(csv_file_path_list)} CSV files in {input_folder}")

    # assign functionality to via smarts 
    amine_pattern = Chem.MolFromSmarts('[NH2]')
    aldehyde_pattern = Chem.MolFromSmarts('[CX3H1](=O)[#6]')

    # iterate through the CSV files in the csv_file_path_list
    for file in csv_file_path_list:
        # extract the base name of the file which doesnt include the extension _peak_* 
        csv_base_name = file.split('_peak_')[0]
        csv_peak_number = file.split('_peak')[1].replace('_refined.csv', '')  # Extract the peak number and remove the .csv extension
        name = csv_base_name + csv_peak_number # Create a name for the output files based on the base name and peak number
    # for each reaction in reactions, extract the relevant data
        for reaction in reactions:
            # Check if the base name matches the reaction's base_name
            reaction_base_name = reaction.get('base_name') # to find which reaction this CSV peak belongs to
            if csv_base_name == reaction_base_name:
                if reaction.get('num_components') == 2:
                    smiles_aldehyde = reaction.get('smiles_aldehyde')
                    topicity_aldehyde = len(Chem.MolFromSmiles(smiles_aldehyde).GetSubstructMatches(aldehyde_pattern)) # number of reactive sites in the aldehyde
                    smiles_amine1 = reaction.get('smiles_amine')
                    topicity_amine1 = len(Chem.MolFromSmiles(smiles_amine1).GetSubstructMatches(amine_pattern)) # number of reactive sites in the amine
                    smiles_amine2 = None
                    topicity_amine2 = 0

                    # as two component, use poc_find_solutions and poc_calc_formulas
                    precursor_combinations = poc_find_solutions(max_precursors, topicity_aldehyde, topicity_amine1)
                    poc_formula_dict = poc_calc_formulas(precursor_combinations, smiles_aldehyde, smiles_amine1, 6)

                elif reaction.get('num_components') == 3:
                    smiles_aldehyde = reaction.get('smiles_aldehyde')
                    topicity_aldehyde = len(Chem.MolFromSmiles(smiles_aldehyde).GetSubstructMatches(aldehyde_pattern)) # number of reactive sites in the aldehyde
                    smiles_amine1 = reaction.get('smiles_amine1')
                    topicity_amine1 = len(Chem.MolFromSmiles(smiles_amine1).GetSubstructMatches(amine_pattern)) # number of reactive sites in the first amine
                    smiles_amine2 = reaction.get('smiles_amine2')
                    topicity_amine2 = len(Chem.MolFromSmiles(smiles_amine2).GetSubstructMatches(amine_pattern)) # number of reactive sites in the second amine

                    # as three component, use poc_find_solutions_ternary and poc_calc_formulas_ternary
                    precursor_combinations = poc_find_solutions_ternary(max_precursors, topicity_aldehyde, topicity_amine1, topicity_amine2)
                    poc_formula_dict = poc_calc_formulas_ternary(precursor_combinations, smiles_aldehyde, smiles_amine1, smiles_amine2, 6)
                else:
                    raise ValueError("Unsupported number of components in reaction. Only 2 or 3 components are supported.")
                
                # Process each reaction to get the isotopes and matching peaks
                isotopes_dict = get_top_10_isotopes(poc_formula_dict)
                mz_intensity_data = read_csv(os.path.join(input_folder, file)) # Read the refined CSV file and extract the m/z and intensity data
                matching_peaks_df = find_matching_isotopes_df(mz_intensity_data, isotopes_dict) # Find matching isotopes in the data
                if matching_peaks_df.empty:
                    print(f"Skipping {name}: No matching peaks found.")
                    continue

                result_df = calculate_mz_differences(matching_peaks_df) # Calculate m/z differences and infer charge from the matching peaks
                result_filtered_df = result_df[result_df['predicted_charge'] == result_df['charge']] # Filter results to keep only those with matching charge
                result_filtered_df = result_filtered_df.sort_values(by='found_intensity', ascending=False) # redorder to highest intensity first

                # Save filtered results
                output_file = os.path.join(output_folder, f"{name}_results.csv")
                result_filtered_df.to_csv(output_file, index=False)
                print(f"Saved: {output_file}")

                # Reload and prepare plot data from the original CSV file
                output_file = os.path.join(input_folder, file)
                df = pd.read_csv(output_file)
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
                smiles_list = [smiles_aldehyde, smiles_amine1] + ([smiles_amine2] if smiles_amine2 else [])
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
                ax.set_title(f"Matching Peaks: {name}")
                plt.tight_layout()
                plt.savefig(os.path.join(output_folder, f"{name}_plot.png"), dpi=300)
                plt.close()

                # Save processed data
                df2_file =  os.path.join(output_folder, f"{name}_processed.csv")
                df2.to_csv(df2_file, index=False)
            else:
                continue  # Skip to the next CSV file if no matching reaction is found
    print("Processing complete.")
if __name__ == "__main__":
    main()
    print("All files processed successfully.")
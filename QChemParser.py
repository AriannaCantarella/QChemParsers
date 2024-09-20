import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scienceplots

class QChemParser:
    def __init__(self, file_path):
        self.file_path = file_path
        self.file_lines = self._read_lines()
        self.file_content = self._read_file()
 
    #should differentiate DFT from HF methods, but still under development
        
    def _read_lines(self):
        """Helper method to read the file content."""
        with open(self.file_path, 'r') as file:
            lines = file.readlines()
            return lines
        
    def _read_file(self):
        """Another helper method to read the file content."""
        with open(self.file_path, 'r') as file:
            return file.read()
        
    ### Only HF methods ###

    def parse_soc_output_RASCI(self):
        """saves data in the  Interstate Transition Properties section from RASCI (and possibly others) calculations"""
        data = []
        current_data = {}
        soc_matrix = {}
        soc_matrix_started = False
        soc_columns = []

        for line in self.file_lines:
            # Detect the start of a new state block
            if line.startswith("State A: Root"):
                # Append previous data before resetting
                if current_data:
                    current_data['soc_matrix'] = soc_matrix  # Ensure SOC matrix is saved
                    data.append(current_data)

                # Reset variables for the new block
                current_data = {}
                soc_matrix = {}
                soc_matrix_started = False
                soc_columns = []

                # Extract Root A
                root_a = int(line.split()[-1])
                current_data["root_A"] = root_a
            
            elif line.startswith("State B: Root"):
                # Extract Root B
                root_b = int(line.split()[-1])
                current_data["root_B"] = root_b

            elif line.startswith("KET: S',Sz'"):
                # Extract spin multiplicity and Sz' for Root A
                parts = re.findall(r"Root\s*(\d+)\).*=\s*(\d+\.?\d*)\s*(\d+\.?\d*)", line)
                if parts:
                    current_data["spin_A"], current_data["sz_A"] = float(parts[0][1]), float(parts[0][2])

            elif line.startswith("BRA: S ,Sz"):
                # Extract spin multiplicity and Sz for Root B
                parts = re.findall(r"Root\s*(\d+)\).*=\s*(\d+\.?\d*)\s*(\d+\.?\d*)", line)
                if parts:
                    current_data["spin_B"], current_data["sz_B"] = float(parts[0][1]), float(parts[0][2])

            elif line.startswith("Clebsh-Gordan coefficient"):
                # Handle cases where SOC elements are skipped
                if "Skipping SOCs" in self.file_content[self.file_content.index(line) + 2]:  # Check if SOCs are skipped
                    current_data["soc_elements"] = None
                    current_data['soc_matrix'] = None
                else:
                    current_data["soc_elements"] = True

            elif line.startswith(" 1-elec SOC matrix"):
                # Start SOC matrix extraction
                soc_matrix_started = True
                soc_columns = []  # Reset columns each time we start a new SOC matrix extraction
            
            # Extract column headers (Sz' values) on the first line of SOC matrix
            elif soc_matrix_started and re.search(r"\|Sz'= *([-\d.]+)>", line.strip()):
                soc_columns = re.findall(r"\|Sz'= *([-\d.]+)>", line.strip())

            
            # Extract SOC matrix rows: <Sz=X| ...> values
            elif soc_matrix_started and re.match(r" *<Sz=\s*([-\d.]+)\|", line.strip()):
                match = re.match(r" *<Sz=\s*([-\d.]+)\|(.+)", line.strip())
                if match:
                    sz_row = match.group(1).strip()  # Extract Sz row value
                    elements = match.group(2).split()  # SOC matrix elements
                    # Process each pair of real and imaginary numbers
                    if len(elements) == 2:
                        real_part = float(elements[0])
                        imag_part = float(elements[1].replace('i', ''))
                        sz_col = soc_columns[0]

                            # Initialize the row if not already
                        if sz_col not in soc_matrix:
                            soc_matrix[sz_col] = {}

                        soc_matrix[sz_col][sz_row] = complex(real_part, imag_part)

                    else:
                        soc_values = []
                        soc_matrix[sz_row] = {}
                        for idx in range(0, len(elements), 2):
                            real_part = float(elements[idx])
                            imag_part = float(elements[idx + 1].replace('i', ''))
                            soc_values.append((real_part,imag_part))
                            sz_col = soc_columns[idx // 2]  # Match Sz' value from the column headers
                            soc_matrix[sz_row][sz_col] = complex(real_part, imag_part)




            elif line.startswith("1-elec SOCC ="):
                # Extract the 1-elec SOCC value
                one_elec_socc = float(line.split('=')[-1].strip().split()[0])
                current_data["one_elec_SOCC"] = one_elec_socc

            elif line.startswith("**************************************************"):
                # Handle the end of a state block
                if current_data:
                    current_data['soc_matrix'] = soc_matrix  # Ensure SOC matrix is saved
                    data.append(current_data)
                    current_data = {}
                    soc_matrix = {}
                    soc_matrix_started = False

        # Append the last data block
        if current_data:
            current_data['soc_matrix'] = soc_matrix  # Ensure SOC matrix is saved
            data.append(current_data)

        return data
    
    def parse_soc_output_RASCI_multiplematrices(self):
        """Saves data in the Interstate Transition Properties section from RASCI calculations."""
        data = []
        current_data = {}
        soc_matrix_started = False
        soc_columns = []
        current_matrix_label = None

        for line in self.file_lines:
            # Detect the start of a new state block
            if line.startswith("State A: Root"):
                # Append previous data before resetting
                if current_data:
                    data.append(current_data)

                # Reset variables for the new block
                current_data = {}
                soc_matrix_started = False
                soc_columns = []
                current_matrix_label = None
                current_data["1-elec_SOC_matrix"] = {}
                current_data["2-elec_SOC_matrix"] = {}
                current_data["Total_SOC_matrix"] = {}

                # Extract Root A
                root_a = int(line.split()[-1])
                current_data["root_A"] = root_a
            
            elif line.startswith("State B: Root"):
                # Extract Root B
                root_b = int(line.split()[-1])
                current_data["root_B"] = root_b

            elif line.startswith("KET: S',Sz'"):
                # Extract spin multiplicity and Sz' for Root A
                parts = re.findall(r"Root\s*(\d+)\).*=\s*(\d+\.?\d*)\s*(\d+\.?\d*)", line)
                if parts:
                    current_data["spin_A"], current_data["sz_A"] = float(parts[0][1]), float(parts[0][2])

            elif line.startswith("BRA: S ,Sz"):
                # Extract spin multiplicity and Sz for Root B
                parts = re.findall(r"Root\s*(\d+)\).*=\s*(\d+\.?\d*)\s*(\d+\.?\d*)", line)
                if parts:
                    current_data["spin_B"], current_data["sz_B"] = float(parts[0][1]), float(parts[0][2])

            elif line.startswith("Clebsh-Gordan coefficient"):
                # Handle cases where SOC elements are skipped
                if "Skipping SOCs" in self.file_content[self.file_content.index(line) + 2]:  # Check if SOCs are skipped
                    current_data["soc_elements"] = None
                else:
                    current_data["soc_elements"] = True

            # Detect different SOC matrix types
            elif line.startswith(" 1-elec SOC matrix"):
                # Start SOC matrix extraction for 1-electron SOC matrix
                soc_matrix_started = True
                soc_columns = []  # Reset columns each time we start a new SOC matrix extraction
                current_matrix_label = "1-elec_SOC_matrix"  # Label for 1-electron SOC matrix
                current_data[current_matrix_label] = {}

            elif line.startswith(" 2-elec mean-field SOC matrix"):
                # Start SOC matrix extraction for 2-electron SOC matrix
                soc_matrix_started = True
                soc_columns = []  # Reset columns
                current_matrix_label = "2-elec_SOC_matrix"  # Label for 2-electron SOC matrix
                current_data[current_matrix_label] = {}

            elif line.startswith(" Total mean-field SOC matrix"):
                # Start SOC matrix extraction for total SOC matrix
                soc_matrix_started = True
                soc_columns = []  # Reset columns
                current_matrix_label = "Total_SOC_matrix"  # Label for total SOC matrix
                current_data[current_matrix_label] = {}

            # Extract column headers (Sz' values) on the first line of SOC matrix
            elif soc_matrix_started and re.search(r"\|Sz'= *([-\d.]+)>", line.strip()):
                soc_columns = re.findall(r"\|Sz'= *([-\d.]+)>", line.strip())

            # Extract SOC matrix rows: <Sz=X| ...> values
            elif soc_matrix_started and re.match(r" *<Sz=\s*([-\d.]+)\|", line.strip()):
                match = re.match(r" *<Sz=\s*([-\d.]+)\|(.+)", line.strip())
                if match:
                    sz_row = match.group(1).strip()  # Extract Sz row value
                    elements = match.group(2).split()  # SOC matrix elements

                    # Process each pair of real and imaginary numbers
                    if len(elements) == 2:
                        real_part = float(elements[0])
                        imag_part = float(elements[1].replace('i', ''))
                        sz_col = soc_columns[0]

                        # Initialize the row if not already
                        if sz_col not in current_data[current_matrix_label]:
                            current_data[current_matrix_label][sz_col] = {}

                        current_data[current_matrix_label][sz_col][sz_row] = complex(real_part, imag_part)

                    else:
                        for idx in range(0, len(elements), 2):
                            real_part = float(elements[idx])
                            imag_part = float(elements[idx + 1].replace('i', ''))
                            sz_col = soc_columns[idx // 2]  # Match Sz' value from the column headers
                            
                            # Initialize the row if not already
                            if sz_col not in current_data[current_matrix_label]:
                                current_data[current_matrix_label][sz_col] = {}

                            current_data[current_matrix_label][sz_col][sz_row] = complex(real_part, imag_part)

            elif line.startswith("1-elec SOCC ="):
                # Extract the 1-elec SOCC value
                one_elec_socc = float(line.split('=')[-1].strip().split()[0])
                current_data["one_elec_SOCC"] = one_elec_socc

            elif line.startswith("**************************************************"):
                # Handle the end of a state block
                if current_data:
                    data.append(current_data)
                    current_data = {}
                    soc_matrix_started = False

        # Append the last data block
        if current_data:
            data.append(current_data)

        return data


    
    def extract_soc_info(self, soc_data, root_A=None, root_B=None, only_socc=False, extract_matrix=False):
        """Creates dataframes with the SOC output data and provides the SOC matrices."""
        """If extract_matrix=True the method reutrns two DFs: the usual df object and the df with the SOC matrix"""

        filtered_results = []
        soc_matrices = []
        
        # Iterate through the SOC data
        for entry in soc_data:
            # Check if Root A matches (mandatory)
            if root_A is not None and entry["root_A"] != root_A:
                continue
            
            # Check if Root B matches (optional)
            if root_B is not None and entry["root_B"] != root_B:
                continue
            
            # If only SOCC is requested, extract only that value
            if only_socc:
                socc_value = entry["one_elec_SOCC"] if entry["one_elec_SOCC"] is not None else None
                filtered_results.append({
                    "Root_A": entry["root_A"],
                    "Root_B": entry["root_B"],
                    "Spin_A": entry["spin_A"],
                    "Spin_B": entry["spin_B"],
                    "1-elec SOCC (cm-1)": socc_value
                })
            else:
                # Add the full entry if not just SOCC
                filtered_results.append({
                    "Root_A": entry["root_A"],
                    "Root_B": entry["root_B"],
                    "Spin_A": entry["spin_A"],
                    "Spin_B": entry["spin_B"],
                    "SOC_Matrix": entry["soc_matrix"]
                })
            
            # If extract_matrix is True, save the SOC matrix in a separate list
            if extract_matrix and entry["soc_matrix"] != {}:
                df_soc_matrix = pd.DataFrame.from_dict(entry['soc_matrix'])
                if entry['spin_A'] == 0.00:
                    df_soc_matrix.columns = pd.MultiIndex.from_tuples([('Sz A', '0.00')])
                    df_soc_matrix.index.names = ['SzB'] 
                elif entry['spin_B'] == 0.00:
                    df_soc_matrix.columns = pd.MultiIndex.from_tuples([('Sz B', '0.00')])
                    df_soc_matrix.index.names = ['SzA'] 
                soc_matrices.append(df_soc_matrix)  # Move this inside the if block

            if extract_matrix and entry["soc_matrix"] == {}:
                soc_matrices.append(None)
        
        # Convert the filtered results into a pandas DataFrame for better handling
        df_filtered = pd.DataFrame(filtered_results)
        
        # If extract_matrix is True, return both the filtered DataFrame and the list of SOC matrices
        if extract_matrix:
            return df_filtered, soc_matrices
        
        return df_filtered

    def build_soc_matrix(self, soc_data, root_A=None, root_B=None, only_socc=False, extract_matrix=False):
        """Creates dataframes with the SOC output data and provides the SOC matrices."""
        
        filtered_results = []
        keys_quintuplet = [-2, -1, 0, 1, 2]
        keys_triplets = [-1, 0, 1]

        # Iterate through the SOC data
        for entry in soc_data:
            # Check if Root A matches (mandatory)
            if root_A is not None and entry["root_A"] != root_A:
                continue

            # Check if Root B matches (optional)
            if root_B is not None and entry["root_B"] != root_B:
                continue

            # If only SOCC is requested, extract only that value
            if only_socc:
                socc_value = entry["one_elec_SOCC"] if entry["one_elec_SOCC"] is not None else None
                filtered_results.append({
                    "Root_A": entry["root_A"],
                    "Root_B": entry["root_B"],
                    "Spin_A": entry["spin_A"],
                    "Spin_B": entry["spin_B"],
                    "1-elec SOCC (cm-1)": socc_value
                })
            else:
                # If not just SOCC, process the SOC matrix
                if extract_matrix and entry["soc_elements"] == True :

                    # Handling singlet-triplet cases
                    if entry["spin_B"] == 1.0 and entry["spin_A"] == 0.0:
                        for Sz_B, soc_value in entry["1-elec_SOC_matrix"]['0.00'].items():
                            filtered_results.append({
                                "Root_A": entry["root_A"],
                                "Root_B": entry["root_B"],
                                "Spin_A": int(entry["spin_A"]),
                                "Spin_B": int(entry["spin_B"]),
                                "Sz_A": int(0.0),
                                "Sz_B": int(float(Sz_B)),  # Extract Sz_B or Sz_A value
                                "SOC element (cm-1)": soc_value
                            })
                    if entry["spin_A"] == 1.0 and entry["spin_B"] == 0.0:
                        for Sz_B, soc_value in entry["1-elec_SOC_matrix"].items():
                            for Sz_A, value in soc_value.items():
                                filtered_results.append({
                                    "Root_A": entry["root_A"],
                                    "Root_B": entry["root_B"],
                                    "Spin_A": int(entry["spin_A"]),
                                    "Spin_B": int(entry["spin_B"]),
                                    "Sz_A": int(float(Sz_B)),  # Extract Sz_B or Sz_A value
                                    "Sz_B": int(0.0),
                                    "SOC element (cm-1)": value
                                })

                    # Handling singlet/triplet-quintuplets cases
                    if entry["spin_B"] == 2.0 and (entry["spin_A"] == 1.0 or entry['spin_A'] == 0.0): 
                        if entry["1-elec_SOC_matrix"] != {}:
                            for key1, subdict in entry["1-elec_SOC_matrix"].items():
                                for key2, value in subdict.items():
                                    filtered_results.append({
                                        "Root_A": entry["root_A"],
                                        "Root_B": entry["root_B"],
                                        "Spin_A": int(entry["spin_A"]),
                                        "Spin_B": int(entry["spin_B"]),
                                        "Sz_A": int(float(key1)),  # Extract Sz_B or Sz_A value
                                        "Sz_B": int(float(key2)),
                                        "SOC element (cm-1)": value
                                    })
                        else:
                            for key in keys_quintuplet:
                                filtered_results.append({
                                    "Root_A": entry["root_A"],
                                    "Root_B": entry["root_B"],
                                    "Spin_A": int(entry["spin_A"]),
                                    "Spin_B": int(entry["spin_B"]),
                                    "Sz_A": int(0),  # Extract Sz_B or Sz_A value
                                    "Sz_B": int(key),
                                    "SOC element (cm-1)": None 
                                })
                    if entry["spin_A"] == 2.0 and (entry["spin_B"] == 1.0 or entry["spin_B"] == 0.0): 
                        if entry["1-elec_SOC_matrix"] != {}:
                            for key1, subdict in entry["1-elec_SOC_matrix"].items():
                                for key2, value in subdict.items():
                                    filtered_results.append({
                                        "Root_A": entry["root_A"],
                                        "Root_B": entry["root_B"],
                                        "Spin_A": int(entry["spin_A"]),
                                        "Spin_B": int(entry["spin_B"]),
                                        "Sz_A": int(float(key1)),  # Extract Sz_B or Sz_A value
                                        "Sz_B": int(float(key2)),
                                        "SOC element (cm-1)": value
                                    })
                        else:
                            for key in keys_quintuplet:
                                filtered_results.append({
                                    "Root_A": entry["root_A"],
                                    "Root_B": entry["root_B"],
                                    "Spin_A": int(entry["spin_A"]),
                                    "Spin_B": int(entry["spin_B"]),
                                    "Sz_A": key,  # Extract Sz_B or Sz_A value
                                    "Sz_B": 0,
                                    "SOC element (cm-1)": None
                                })

                    # Handling triplet-triplet cases
                    if entry["spin_B"] == 1.0 and entry["spin_A"] == 1.0: 
                                for key1 in keys_triplets:
                                    for key2 in keys_triplets:
                                        filtered_results.append({
                                            "Root_A": entry["root_A"],
                                            "Root_B": entry["root_B"],
                                            "Spin_A": int(entry["spin_A"]),
                                            "Spin_B": int(entry["spin_B"]),
                                            "Sz_A": int(key1),  # Extract Sz_B or Sz_A value
                                            "Sz_B": int(key2),
                                            "SOC element (cm-1)": None
                                        })

                    # Handling singlet-singlet                    
                    if entry["spin_A"] == 0.0 and entry["spin_B"] == 0.0: 
                        filtered_results.append({
                            "Root_A": entry["root_A"],
                            "Root_B": entry["root_B"],
                            "Spin_A": int(entry["spin_A"]),
                            "Spin_B": int(entry["spin_B"]),
                            "Sz_A": 0,  # Extract Sz_B or Sz_A value
                            "Sz_B": 0,
                            "SOC element (cm-1)": None
                        })
                else:
                    # If there's no SOC matrix, we can still add the basic info
                    filtered_results.append({
                        "Root_A": entry["root_A"],
                        "Root_B": entry["root_B"],
                        "Spin_A": int(entry["spin_A"]),
                        "Spin_B": int(entry["spin_B"]),
                        "Sz_A": int(float(entry["sz_A"])),
                        "Sz_B": int(float(entry['sz_B'])),
                        "SOC element (cm-1)": None
                    })

        # Convert the filtered results into a pandas DataFrame for better handling
        df_filtered = pd.DataFrame(filtered_results)

        return df_filtered


    def rasci_amplitude_coefficients(file_path):
        """
        Processes a RAS-CI output file to extract the amplitude matrix and excitation energies.

        Parameters:
        file_path (str): Path to the RAS-CI output file.

        Returns:
        amplitude_df (pd.DataFrame): Transposed DataFrame containing the amplitude matrix, with states as rows and 'alpha|beta' as columns.
        excitation_energies (list): List containing excitation energies corresponding to each state.
        """
        amplitude_df = pd.DataFrame()
        excitation_energies = []
        alpha_beta_order = []  # List to keep track of alpha|beta order

        with open(file_path, 'r') as f:
            lines = f.readlines()

        in_state_block = False
        current_state = None
        current_energy = None

        for i, line in enumerate(lines):
            if "RAS-CI total energy for state" in line:
                current_state = int(re.search(r"state\s+(\d+):", line).group(1))
                in_state_block = True

            if in_state_block and "Excitation energy (eV) =" in line:
                current_energy = float(re.search(r"Excitation energy \(eV\) =\s+([\d\.\-E]+)", line).group(1))
                excitation_energies.append(current_energy)

            if in_state_block and "AMPLITUDE" in line:
                start_index = i + 2
                for j in range(start_index, len(lines)):
                    if "--------------------------------------------------" in lines[j]:
                        end_index = j
                        break
                
                # Extract amplitude data
                current_amplitudes = []
                for k in range(start_index, end_index):
                    if "|" in lines[k] and "AMPLITUDE" not in lines[k]:
                        parts = [p.strip() for p in lines[k].split("|") if p.strip()]
                        if len(parts) == 3:
                            alpha = parts[0]
                            beta = parts[1]
                            try:
                                amplitude = float(parts[2])
                                alpha_beta = f"{alpha}{beta}"
                                current_amplitudes.append((alpha_beta, amplitude))
                                if alpha_beta not in alpha_beta_order:
                                    alpha_beta_order.append(alpha_beta)  # Store order
                            except ValueError:
                                continue
                
                # Create DataFrame for current state and merge
                if current_amplitudes:
                    df_state = pd.DataFrame(current_amplitudes, columns=['alpha|beta', f'State {current_state}'])
                    df_state.set_index('alpha|beta', inplace=True)
                    amplitude_df = pd.merge(amplitude_df, df_state, left_index=True, right_index=True, how='outer')

                in_state_block = False  # Reset after processing

        amplitude_df.fillna(0.0, inplace=True)

        # Reorder the columns in the DataFrame based on the order in alpha_beta_order
        amplitude_df = amplitude_df.T[alpha_beta_order]

        return amplitude_df, excitation_energies


    def reading_MOs(self, unrestricted=False, n_th_occurrence=-1):
        """Reads and parses the last occurence of orbital sets. Pay attention to when many calculations are printed in the same output file
           Juat added the possibility to select the n-th occurence in the output file. """
        # Read the content of the file

        content = self.file_content

        au_to_ev = 27.211324570273

        # Find all occurrences of the start and end indices of the sections
        start_indices = [match.start() for match in re.finditer(r"Orbital Energies \(a\.u\.\)", content)]
        end_indices = [match.end() for match in re.finditer(r"Ground-State Mulliken Net Atomic Charges", content)]

        if not start_indices or not end_indices:
            raise ValueError("Orbital Energies or Mulliken Charges section not found in the file.")
     
        if start_indices and end_indices:
            # Select the last occurrence
            start_index = start_indices[n_th_occurrence]
            end_index = end_indices[n_th_occurrence]
            section_content = content[start_index:end_index]

            # Check if both start and end indices are found
            if start_index != -1 and end_index != -1:
                # Extract the section between start and end indices
                orbital_energies_section = content[start_index:end_index]

                if unrestricted:
                    # Unrestricted case: Separate alpha and beta MOs
                    alpha_section = re.search(r"Alpha MOs(.+?)Beta MOs", orbital_energies_section, re.DOTALL)
                    beta_section = re.search(r"Beta MOs(.+?)(Ground-State Mulliken Net Atomic Charges|$)", orbital_energies_section, re.DOTALL)

                    if alpha_section and beta_section:
                        alpha_content = alpha_section.group(1)
                        beta_content = beta_section.group(1)

                        # Extract occupied and virtual values for alpha orbitals
                        alpha_occupied_line = re.search(r"-- Occupied --(.+?)(-- Virtual --|$)", alpha_content, re.DOTALL)
                        alpha_virtual_line = re.search(r"-- Virtual --(.+?)$", alpha_content, re.DOTALL)

                        if alpha_occupied_line and alpha_virtual_line:
                            alpha_occupied_values = [float(match) * au_to_ev for match in re.findall(r"[-+]?\d*\.\d+|\d+", alpha_occupied_line.group(1))]
                            alpha_virtual_values = [float(match) * au_to_ev for match in re.findall(r"[-+]?\d*\.\d+|\d+", alpha_virtual_line.group(1))]
                        else:
                            print("Alpha occupied or virtual section not found.")
                            return None

                        # Extract occupied and virtual values for beta orbitals
                        beta_occupied_line = re.search(r"-- Occupied --(.+?)(-- Virtual --|$)", beta_content, re.DOTALL)
                        beta_virtual_line = re.search(r"-- Virtual --(.+?)$", beta_content, re.DOTALL)

                        if beta_occupied_line and beta_virtual_line:
                            beta_occupied_values = [float(match) * au_to_ev for match in re.findall(r"[-+]?\d*\.\d+|\d+", beta_occupied_line.group(1))]
                            beta_virtual_values = [float(match) * au_to_ev for match in re.findall(r"[-+]?\d*\.\d+|\d+", beta_virtual_line.group(1))]
                        else:
                            print("Beta occupied or virtual section not found.")
                            return None

                        alpha_gap = alpha_virtual_values[0] - alpha_occupied_values[-1]
                        beta_gap = beta_virtual_values[0] - beta_occupied_values[-1]

                        print(f'Alpha gap is {alpha_gap} eV')
                        print(f'Beta gap is {beta_gap} eV')

                        return alpha_occupied_values, alpha_virtual_values, alpha_gap, beta_occupied_values, beta_virtual_values, beta_gap
                    else:
                        print("Alpha or Beta section not found.")
                        return None
                else:
                    # Restricted case
                    # Use regular expressions to find all numerical values in the extracted section
                    values_pattern = r"[-+]?\d*\.\d+|\d+"
                    values = [float(match) for match in re.findall(values_pattern, orbital_energies_section)]

                    # Find the lines containing "Occupied" and "Virtual"
                    occupied_line = re.search(r"-- Occupied --(.+?)--", orbital_energies_section, re.DOTALL)
                    virtual_line = re.search(r"-- Virtual --(.+?)--", orbital_energies_section, re.DOTALL)

                    if occupied_line and virtual_line:
                        # Extract the values from the matched lines
                        occupied_values = [float(match) * au_to_ev for match in re.findall(values_pattern, occupied_line.group(1))]
                        virtual_values = [float(match) * au_to_ev for match in re.findall(values_pattern, virtual_line.group(1))]

                        gap = virtual_values[0] - occupied_values[-1]
                        print(f'Gap is {gap} eV')

                        return occupied_values, virtual_values, gap
                    else:
                        print("Occupied or Virtual section not found.")
                        return None
        else:
            # Handle case where no sections are found
            print("No valid sections found in the file.")
            return None


    
    ### Only DFT methods ###

    def meaningful_quantities_TDDFT(self):

        """ function for extracting excitation_energies, states, multiplicities, 
            strengths, transitions, exc_state_energy from a TDDFT calc output by QChem """
        
    # inputs: - qchem.out file
    #         - the KS orbital energies given by the same calculation

        output_text = self.file_content

        def eV_to_lambda(energy):
            planck_h = 6.6261e-34 #J·s
            c = 299792458 #m/s
            wavelength = planck_h*c/(energy*1.602176565e-19) * 1e9
            return wavelength
        
        # Define regex patterns to extract the desired information
        pattern_state = re.compile(r'Excited state\s+(\d+): excitation energy \(eV\) =\s+([\d.]+)')
        pattern_multiplicity = re.compile(r'Multiplicity:\s+(Singlet|Triplet)')
        pattern_strength = re.compile(r'Strength\s+:\s+([\d.]+)')
        pattern_transition = re.compile(r'D\(\s*(\d+)\)\s+-->\s+V\(\s*(\d+)\)\s+amplitude\s+=\s+(-?\d+\.\d+)')

        # Initialize lists to store the extracted information
        excitation_energies = []
        states = []
        multiplicities = []
        strengths = []
        transitions = []
        states_multipl = []

        # Iterate over each match in the text
        for match in re.finditer(pattern_state, output_text):
            # Extract the state number and excitation energy from the matched groups
            state_number = int(match.group(1))
            energy = float(match.group(2))
            excitation_energies.append(energy)
            states.append(state_number)

        # Iterate over each match in the text
        i = 0  # Counter for Triplets
        j = 0  # Counter for Singlets

        # Iterate over each match in the text
        for match in re.finditer(pattern_multiplicity, output_text):
            # Extract the multiplicity from the matched group
            multiplicity = match.group(1)
            multiplicities.append(multiplicity)
            if multiplicity == 'Triplet':
                i += 1  # Increment the Triplet counter
                states_multipl.append('T' + str(i))
            elif multiplicity == 'Singlet':  # Using elif to make it mutually exclusive
                j += 1  # Increment the Singlet counter
                states_multipl.append('S' + str(j))


        # Iterate over each match in the text
        for match in re.finditer(pattern_strength, output_text):
            # Extract the strength from the matched group
            strength = float(match.group(1))
            strengths.append(strength)

        # Iterate over each match in the text
        for match in re.finditer(pattern_transition, output_text):
            # Extract the transition details from the matched groups
            start_state = int(match.group(1))
            end_state = int(match.group(2))
            amplitude = float(match.group(3))
            transitions.append((start_state, end_state, amplitude))

        # Group transitions by state number
        transition_dict = {}
        for state, transition in zip(states, transitions):
            if state not in transition_dict:
                transition_dict[state] = []
            transition_dict[state].append(transition)


        # Create a DataFrame from the extracted information
        data = []
        for state in states:
            data.append({
                'State': states_multipl[state-1],
                'Excitation Energy (eV)': excitation_energies[state-1],
                'Wavelength' : eV_to_lambda(excitation_energies[state-1]),
                'Multiplicity': multiplicities[state-1],
                'Strength': strengths[state-1],
                'Transitions': transition_dict.get(state, [])
            })
        df = pd.DataFrame(data)
        df.index = pd.RangeIndex(start=1, stop=(len(df) + 1), step=1)


        #df = df[df['Wavelength'] > 300]
        #df['Strength'] = df['Strength'] / max(df['Strength']) * 0.8 #normalizzando dopo il filtering delle wavelength
        # Display the DataFrame
        #print(df)

        return df
    
#ex 
# file_path = '/Users/acantarella/CISS/RAS-CI/test_SOC_parser.out'

# from QChemParser import QChemParser

# parser = QChemParser(file_path)
# soc_output = parser.parse_soc_output_RASCI()
# a, b = parser.extract_soc_info(soc_output, root_A=1, extract_matrix=True)
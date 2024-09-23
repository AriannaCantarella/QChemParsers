# How to use the parser for different tools:

## parsing of the soc matrix from RAS-CI calculation ##

from QChem_Parsers.QChemParser import QChemParser 

filepath = 'your_file_path'
parser = QChemParser(filepath)

soc_data = parser.parse_soc_output_RASCI_multiplematrices() #it builds a dictionary with info about all the states and the soc_matrix

df = parser.build_soc_matrix(soc_data, root_A=a, extract_matrix=True) #one can specify both root_A and root_B

## parsing of the coefficients (amplitudes) in RAS-CI ##

parser = QChemParser(filepath)

coeff_configs, en = parser.rasci_amplitude_coefficients()


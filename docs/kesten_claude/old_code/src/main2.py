import json
import sys

def parse_nasa_thermo(file_content):
    """
    Parses a NASA thermochemistry file (like the attached nasa9.dat).

    This parser is built for the format observed in the file:
    - A 'thermo' keyword.
    - Species blocks, each starting with a species name line.
    - An info line specifying the number of temperature ranges.
    - For each range:
        - 1 line for T-range, 8 exponents, and H(0).
        - 2 lines for 9 coefficients.
    - An 'END' keyword.
    """
    species_data = {}
    lines = file_content.splitlines()
    lines_iter = iter(lines)
    in_thermo_section = False

    current_species_data = None
    species_name = None
    line_number = 0

    try:
        while True:
            line_number += 1
            # Read a line and scrub the non-breaking space (u00a0)
            line = next(lines_iter).replace('\u00a0', ' ')

            # Skip comments
            if line.strip().startswith('!'):
                continue
            
            # Start of data
            if 'thermo' in line.lower():
                in_thermo_section = True
                # Skip the global temp line that follows 'thermo'
                line_number += 1
                next(lines_iter) 
                continue
            
            # End of file
            if 'end' in line.lower():
                break
                
            if not in_thermo_section:
                continue

            # --- Start of a new species entry ---

            # Line 1: Species Name and Comment
            species_name = line[0:18].strip()
            comment = line[18:80].strip()
            
            current_species_data = {
                'species': species_name,
                'comment': comment,
                'ranges': []
            }

            # Line 2: Info Line
            line_number += 1
            info_line = next(lines_iter).replace('\u00a0', ' ')
            
            # Extract number of temperature ranges
            num_ranges = int(info_line[0:2].strip())
            
            # Extract date
            current_species_data['date'] = info_line[3:9].strip()

            # Extract elemental composition (up to 5 elements)
            atoms = {}
            for i in range(5):
                start = 10 + (i * 5)
                atom_info = info_line[start:start+5].strip()
                if atom_info:
                    # Find where the element symbol ends and count begins
                    elem = ""
                    count = ""
                    for char in atom_info:
                        if char.isalpha():
                            elem += char
                        else:
                            count += char
                    if elem and count:
                         atoms[elem] = float(count)
            current_species_data['atoms'] = atoms

            # Extract phase (0=gas, 1=condensed)
            current_species_data['phase'] = int(info_line[30:31].strip())
            
            # Extract molecular weight
            current_species_data['molecular_weight'] = float(info_line[55:71].strip().replace('D', 'E'))
            
            # Extract heat of formation at 298.15K (if available on this line)
            try:
                current_species_data['hf_298_j_mol'] = float(info_line[71:80].strip().replace('D', 'E'))
            except ValueError:
                current_species_data['hf_298_j_mol'] = None

            # Loop over the temperature ranges for this species
            for _ in range(num_ranges):
                
                # T-Range Line (Line 3 of range)
                line_number += 1
                t_line = next(lines_iter).replace('\u00a0', ' ')
                
                t_low = float(t_line[0:12].strip())
                t_high = float(t_line[12:24].strip())
                num_coeffs = int(float(t_line[24:30].strip())) # Should be 7 for this file
                
                # *** FIX 2: Loop 8 times (30 to 70) for 8 exponents ***
                # Exponents (8 terms)
                exponents = []
                for i in range(30, 70, 5): # 30, 35, 40, 45, 50, 55, 60, 65
                    exponents.append(float(t_line[i:i+5].strip()))

                try:
                    # H(0) - H(298.15) in J/mol (or similar reference)
                    h_ref = float(t_line[75:88].strip().replace('D', 'E'))
                except (ValueError, IndexError):
                    h_ref = 0.0

                range_data = {
                    't_low': t_low,
                    't_high': t_high,
                    'num_coeffs': num_coeffs,
                    'exponents': exponents,
                    'h_ref': h_ref
                }

                # Coefficient Lines (must replace 'D' with 'E' for Python's float)
                line_number += 1
                coeff_line_1 = next(lines_iter).replace('\u00a0', ' ').replace('D', 'E')
                line_number += 1
                coeff_line_2 = next(lines_iter).replace('\u00a0', ' ').replace('D', 'E')

                # Extract coefficients using fixed 15-character widths
                coeffs1 = [float(coeff_line_1[i:i+15]) for i in range(0, 75, 15)]
                coeffs2 = [float(coeff_line_2[i:i+15]) for i in range(0, 60, 15)]
                
                all_coeffs = coeffs1 + coeffs2
                
                # This file format provides 7 polynomial coefficients (a1-a7)
                # and 2 integration constants (b1, b2)
                if len(all_coeffs) == 9:
                    range_data['coeffs_a'] = all_coeffs[0:7]
                    range_data['coeffs_b'] = all_coeffs[7:9]
                else:
                    print(f"Warning: Unexpected coefficient count for {species_name}")

                current_species_data['ranges'].append(range_data)
                
            species_data[species_name] = current_species_data

    except StopIteration:
        # Expected end of file
        pass
    except Exception as e:
        # Handle any parsing errors
        print(f"--- ERROR ---")
        print(f"Error parsing species: {species_name} (near line {line_number})")
        print(f"Last species line read: '{line}'")
        print(f"Error message: {e}")
        print(f"---------------")
        pass

    return species_data

if __name__ == "__main__":
    file_name = 'nasa9.dat'
    try:
        # Specify encoding as 'utf-8' to handle special characters
        with open(file_name, 'r', encoding='utf-8') as f:
            file_content = f.read()
        
        print(f"Parsing {file_name}...")
        thermo_data = parse_nasa_thermo(file_content)
        print(f"Successfully parsed {len(thermo_data)} species.")
        
        # --- Example Usage ---
        
        # 1. Print the first 3 parsed species as a JSON object
        print("\n--- Example: First 3 Species Data ---")
        count = 0
        for species_name, data in thermo_data.items():
            if count < 3:
                print(json.dumps(data, indent=2))
                count += 1
            else:
                break
        
        # 2. Get specific data for one species (e.g., Ag)
        if 'Ag' in thermo_data:
            print("\n--- Example: Data for 'Ag' ---")
            print(f"Species: {thermo_data['Ag']['species']}")
            print(f"MW: {thermo_data['Ag']['molecular_weight']}")
            print(f"Phase: {'Gas' if thermo_data['Ag']['phase'] == 0 else 'Condensed'}")
            print(f"Number of T-Ranges: {len(thermo_data['Ag']['ranges'])}")
            
            # Print exponents for the second temperature range of Ag
            ag_range_2 = thermo_data['Ag']['ranges'][1]
            print("\nExponents for 'Ag' (Range 2: 1000-6000K):")
            print(f"  Exponents: {ag_range_2['exponents']}")
            
        # 3. Get specific data for H2O
        if 'H2O' in thermo_data:
             print("\n--- Example: Data for 'H2O' ---")
             print(f"Species: {thermo_data['H2O']['species']}")
             print(f"MW: {thermo_data['H2O']['molecular_weight']}")
             # Get coefficients for the first range (low temp)
             h2o_low_t = thermo_data['H2O']['ranges'][0]
             print(f"Low-T Range ({h2o_low_t['t_low']}K - {h2o_low_t['t_high']}K):")
             print(f"  Exponents: {h2o_low_t['exponents']}")


    except FileNotFoundError:
        print(f"Error: The file '{file_name}' was not found.")
        print("Please make sure the script is in the same directory as the data file.")
    except Exception as e:
        print(f"An error occurred: {e}")
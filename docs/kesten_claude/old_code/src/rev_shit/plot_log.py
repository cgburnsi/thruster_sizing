# plot_log.py
import matplotlib.pyplot as plt
import re
import numpy as np

def parse_log_file(log_path="simulation.log"):
    """
    Parses the full simulation.log file and extracts all data.
    """
    print(f"Parsing '{log_path}'...")
    
    # We'll store data for each region
    regions = {
        'liquid': {'z': [], 'temp': [], 'h': []},
        'lqvp': {'z': [], 'temp': [], 'h': [], 'wfv': []},
        'vapor': {'z': [], 'temp': [], 'h': [], 'c1': [], 'c2': [], 'c3': [], 'c4': [], 'frac3d': []}
    }

    # Regex to match the three main data lines in your log
    re_liq = re.compile(r"Z:(\S+) \[ft\], TEMP:(\S+) \[degR\], H:(\S+) \[BTU/lb\], DHDZ:(\S+) \[ft\]")
    re_lqvp = re.compile(r"Z:(\S+) \[ft\], TEMP:(\S+) \[degR\], H:(\S+) \[BTU/lb\], WFV:\s*(\S+)")
    re_vap = re.compile(r"Z:(\S+) \[ft\], T:(\S+) \[R\], P:(\S+) \[psia\], H:(\S+) \[BTU/lb\], C4:(\S+), C3:(\S+), C2:(\S+), C1:(\S+), FRAC3D:(\S+)")

    with open(log_path, 'r') as f:
        for line in f:
            match_liq = re_liq.match(line)
            match_lqvp = re_lqvp.match(line)
            match_vap = re_vap.match(line)

            if match_liq:
                data = [float(v) for v in match_liq.groups()]
                regions['liquid']['z'].append(data[0])
                regions['liquid']['temp'].append(data[1])
                regions['liquid']['h'].append(data[2])
            
            elif match_lqvp:
                data = [float(v) for v in match_lqvp.groups()]
                regions['lqvp']['z'].append(data[0])
                regions['lqvp']['temp'].append(data[1])
                regions['lqvp']['h'].append(data[2])
                regions['lqvp']['wfv'].append(data[3])

            elif match_vap:
                data = [float(v) for v in match_vap.groups()]
                regions['vapor']['z'].append(data[0])
                regions['vapor']['temp'].append(data[1])
                regions['vapor']['h'].append(data[3]) # group 2 is P
                regions['vapor']['c4'].append(data[4])
                regions['vapor']['c3'].append(data[5])
                regions['vapor']['c2'].append(data[6])
                regions['vapor']['c1'].append(data[7])
                regions['vapor']['frac3d'].append(data[8])
                
    print("Parsing complete.")
    return regions

def plot_full_profile(data):
    """
    Generates three summary plots from the parsed log data.
    """
    
    # --- Plot 1: Full Temperature Profile ---
    
    # Combine data for a continuous plot
    z_full = data['liquid']['z'] + data['lqvp']['z'] + data['vapor']['z']
    t_full = data['liquid']['temp'] + data['lqvp']['temp'] + data['vapor']['temp']

    plt.figure(1, figsize=(10, 6))
    plt.plot(z_full, t_full, 'b-')
    plt.title('Axial Temperature Profile', fontsize=16)
    plt.xlabel('Axial Position (Z) [ft]')
    plt.ylabel('Temperature [degR]')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    
    # Add vertical lines to show the regions
    z_boil_start = data['liquid']['z'][-1]
    z_boil_end = data['lqvp']['z'][-1]
    plt.axvline(x=z_boil_start, color='k', linestyle=':', linewidth=2, label=f'Boiling Starts (Z={z_boil_start:.4f} ft)')
    plt.axvline(x=z_boil_end, color='k', linestyle='--', linewidth=2, label=f'Vapor Region Starts (Z={z_boil_end:.4f} ft)')
    
    # Find and plot the peak temperature
    t_peak = max(data['vapor']['temp'])
    z_peak = data['vapor']['z'][np.argmax(data['vapor']['temp'])]
    plt.plot(z_peak, t_peak, 'ro', label=f'Peak Temp ({t_peak:.0f} R)')
    
    plt.legend()
    plt.tight_layout()

    # --- Plot 2: Vapor Phase Mass Concentrations ---
    plt.figure(2, figsize=(10, 6))
    plt.plot(data['vapor']['z'], data['vapor']['c4'], label='N2H4 (C4)', linewidth=2)
    plt.plot(data['vapor']['z'], data['vapor']['c3'], label='NH3 (C3)', linewidth=2)
    plt.plot(data['vapor']['z'], data['vapor']['c2'], label='N2 (C2)', linewidth=2)
    plt.plot(data['vapor']['z'], data['vapor']['c1'], label='H2 (C1)', linewidth=2)
    plt.title('Vapor Phase Mass Concentrations', fontsize=16)
    plt.xlabel('Axial Position (Z) [ft]')
    plt.ylabel('Mass Concentration (C_i) [lb/ft^3]')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()

    # --- Plot 3: Ammonia Dissociation ---
    plt.figure(3, figsize=(10, 6))
    plt.plot(data['vapor']['z'], data['vapor']['frac3d'], 'r-')
    plt.title('Fractional Ammonia Dissociation (FRAC3D)', fontsize=16)
    plt.xlabel('Axial Position (Z) [ft]')
    plt.ylabel('Fractional Dissociation [-]')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()

    print("Displaying plots...")
    plt.show()

if __name__ == "__main__":
    try:
        log_data = parse_log_file("simulation.log")
        plot_full_profile(log_data)
    except Exception as e:
        print(f"\nAn error occurred:")
        print(f"  {e}")
        print("\nMake sure 'simulation.log' is in the same directory.")
        print("If the log is incomplete, parsing may fail.")
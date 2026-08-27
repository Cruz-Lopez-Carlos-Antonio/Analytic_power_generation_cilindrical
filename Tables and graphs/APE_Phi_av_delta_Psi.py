import mpmath
import os

# Set mpmath precision to 50 digits to prevent catastrophic cancellation
mpmath.mp.dps = 50

def calculate_ape():
    """
    Calculates the Absolute Percentage Error (APE) between a semi-analytic reference 
    model and a direct numerical solution.
    """
    file_analytic = "PhiAv_vs_delta_varying_PsiS_analytic.txt"
    file_direct = "Direct_PhiAv_vs_delta_varying_PsiS.txt"
    file_output = "APE_PhiAv_vs_delta_varying_PsiS.txt"

    try:
        with open(file_analytic, 'r') as fa, open(file_direct, 'r') as fd:
            lines_a = fa.readlines()
            lines_d = fd.readlines()
    except FileNotFoundError:
        print("Error: Input files not found.")
        print(f"Ensure that '{file_analytic}' and '{file_direct}' are located in the current directory.")
        return

    # Retrieve headers from the reference (analytic) file
    headers = lines_a[0].strip().split('\t')

    with open(file_output, 'w') as f_out:
        # Write headers preserving the tab-separated format
        f_out.write('\t'.join(headers) + '\n')

        # Iterate through the data rows (skipping the header line)
        for la, ld in zip(lines_a[1:], lines_d[1:]):
            parts_a = la.strip().split('\t')
            parts_d = ld.strip().split('\t')

            # The first element is the 'delta' parameter
            delta = parts_a[0]
            ape_row = [delta]

            # Iterate through the numerical columns
            for val_a_str, val_d_str in zip(parts_a[1:], parts_d[1:]):
                # Convert to mpmath.mpf for multi-precision processing
                val_a = mpmath.mpf(val_a_str)
                val_d = mpmath.mpf(val_d_str)

                # Calculate APE: |(Numerical - Analytic) / Analytic| * 100
                if val_a != 0:
                    ape = mpmath.fabs((val_d - val_a) / val_a) * mpmath.mpf('100')
                else:
                    ape = mpmath.mpf('0')
                
                # Format APE to 5 significant digits using scientific notation (4 decimals)
                ape_formatted = "{:.4e}".format(float(ape))
                ape_row.append(ape_formatted)

            f_out.write('\t'.join(ape_row) + '\n')
            
    print(f"Execution successful. Error table generated at: {file_output}")

if __name__ == '__main__':
    calculate_ape()

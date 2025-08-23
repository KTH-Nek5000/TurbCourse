import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 10,
    'font.family': 'serif',
    'font.serif': 'Times New Roman',
    'mathtext.fontset': 'cm',
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'axes.grid': True,
    'axes.edgecolor': 'black',
    'axes.linewidth': 0.8,
    'grid.linestyle': '--',
    'grid.alpha': 0.7,
    'grid.color': 'gray',
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.size': 5,
    'ytick.major.size': 5,
    'lines.linewidth': 1.5,
    'lines.markersize': 5,
    'legend.fontsize': 10,
    'legend.frameon': False,
    'legend.loc': 'best',
    'savefig.dpi': 300,
    'savefig.format': 'png',
})


def freqPOD(A2,maxModes,nplt,timestep,info):

    #A2 -- POD Coefficients (m,n) where m is the number of snapshots and n is the mode
    #Mode 0 is the time average.

    fileName = info['outputPath']+'PODFreq.txt'
    F = open(fileName,'w')
    F.write("# Snapshots POD result \n") 
    F.write("# Dominant frequencies calculated using FFT \n") 
    F.write("# According to our convention the mode 0 is the mean value \n") 
    F.write("# ------------------------------------------------------------------\n") 
    F.write("# Mode\t F1\t  F2\t  F3\t F4\t F5\n") 
    F.write("# ------------------------------------------------------------------\n") 
    
    # Define frequencies for FFT
    frequencies = np.fft.rfftfreq(nplt, d=timestep)  # Generates only positive frequencies
    print(frequencies)

    for i in range(1,maxModes):
        # Step 1: Compute FFT
        fft_result = np.fft.fft(A2[:, i] - np.mean(A2[:, i]))
        power = np.abs(fft_result[:nplt // 2 + 1])  # Compute power spectrum (positive frequencies only)

        # Step 2: Filter the power spectrum
        mean_power = np.mean(power)
        std_power = np.std(power)
        threshold = mean_power + 2 * std_power
        filtered_power = np.where(power > threshold, power, 0)

        # Step 3: Identify dominant frequencies
        idx = np.argsort(filtered_power)  # Sort indices based on filtered power
        f = frequencies[np.flip(idx)[0:5]]  # Top 5 frequencies
        F.write("%g\t%g\t%g\t%g\t%g\t%g \n" % \
                (i, f[0],f[1],f[2],f[3],f[4])) 

        plt.figure(figsize=(4,3))
        plt.plot(frequencies, power, label='FFT Power', alpha=0.6)
        plt.plot(frequencies, filtered_power, label='Filtered Power', linestyle="solid")
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Power')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(info['outputPath']+f"Mode{i}.png")
        plt.close() 
    
    F.close()

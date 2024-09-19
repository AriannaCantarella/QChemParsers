import matplotlib.pyplot as plt
import numpy as np
import scienceplots

def process_energies(energies, state_numbers, degeneracy_threshold=0.02):
    unique_energies = []
    degenerate_counts = []
    state_numbers_grouped = []

    energies_sorted = sorted(energies)
    i = 0
    while i < len(energies_sorted):
        count = 1
        energy = energies_sorted[i]
        states = [state_numbers[energies.index(energy)]]
        while (i + 1 < len(energies_sorted) and 
               abs(energies_sorted[i + 1] - energy) < degeneracy_threshold):
            i += 1
            count += 1
            states.append(state_numbers[energies.index(energies_sorted[i])])
        unique_energies.append(energy)
        degenerate_counts.append(count)
        state_numbers_grouped.append(states)
        i += 1

    return unique_energies, degenerate_counts, state_numbers_grouped

def plot_energy_spectrum(ax, occupied_values, virtual_values, title, split_alpha_beta=False, 
                         alpha_occupied=None, beta_occupied=None, alpha_virtual=None, beta_virtual=None,
                         degeneracy_threshold=0.02):
    ax.set_title(title)
    ax.yaxis.grid(True, which='both')
    ax.set_ylim([-8, 0])
    ax.xaxis.grid(False)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.yaxis.set_major_locator(plt.MultipleLocator(1))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.2))
    ax.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')

    if split_alpha_beta and alpha_occupied and beta_occupied and alpha_virtual and beta_virtual:
        # Process alpha orbitals
        alpha_state_numbers = [i + 1 for i in range(len(alpha_occupied))]
        alpha_unique_energies, alpha_degenerate_counts, alpha_state_numbers_grouped = process_energies(alpha_occupied, alpha_state_numbers, degeneracy_threshold)

        # Process beta orbitals
        beta_state_numbers = [i + 1 for i in range(len(beta_occupied))]
        beta_unique_energies, beta_degenerate_counts, beta_state_numbers_grouped = process_energies(beta_occupied, beta_state_numbers, degeneracy_threshold)

        # Plotting alpha orbitals
        for energy, count, states in zip(alpha_unique_energies, alpha_degenerate_counts, alpha_state_numbers_grouped):
            if count == 1:
                ax.hlines(energy, xmin=0, xmax=0.48, color='blue', linewidth=1.5)
                ax.text(0.11, energy + 0.075, f'{states[0]}', va='center', clip_on=True, fontsize=8)
            else:
                ax.hlines(energy, xmin=0, xmax=0.21, color='blue', linewidth=1.5)
                ax.hlines(energy, xmin=0.27, xmax=0.48, color='blue', linewidth=1.5)
                for j, state in enumerate(states):
                    offset = (0.1 if j % 2 == 0 else 0.37)
                    ax.text(offset, energy + 0.075, f'{state}', va='center', ha='center', clip_on=True, fontsize=8)

        # Plotting beta orbitals
        for energy, count, states in zip(beta_unique_energies, beta_degenerate_counts, beta_state_numbers_grouped):
            if count == 1:
                ax.hlines(energy, xmin=0.52, xmax=1, color='green', linewidth=1.5)
                ax.text(0.63, energy + 0.075, f'{states[0]}', va='center', clip_on=True, fontsize=8)
            else:
                ax.hlines(energy, xmin=0.52, xmax=0.73, color='green', linewidth=1.5)
                ax.hlines(energy, xmin=0.77, xmax=1, color='green', linewidth=1.5)
                for j, state in enumerate(states):
                    offset = (0.63 if j % 2 == 0 else 0.88)
                    ax.text(offset, energy + 0.075, f'{state}', va='center', ha='center', clip_on=True, fontsize=8)

        # Plotting virtual alpha orbitals
        for energy, num_state in zip(alpha_virtual, [i + 1 for i in range(len(alpha_virtual))]):
            ax.hlines(energy, alpha=1, linewidth=1.5, xmin=0, xmax=0.48, color='firebrick')
            ax.text(0.24, energy + 0.075, f'{num_state}', va='center', clip_on=True, fontsize=8)

        # Plotting virtual beta orbitals
        for energy, num_state in zip(beta_virtual, [i + 1 for i in range(len(beta_virtual))]):
            ax.hlines(energy, alpha=1, linewidth=1.5, xmin=0.52, xmax=1, color='darkorange')
            ax.text(0.76, energy + 0.075, f'{num_state}', va='center', clip_on=True, fontsize=8)
        
        ax.set_xlabel('alpha orbitals      beta orbitals', fontsize=14)

    else:
        # Process the energies to find unique levels and degeneracies
        num_state_occ = [i + 1 for i in range(len(occupied_values))]
        unique_energies, degenerate_counts, state_numbers = process_energies(occupied_values, num_state_occ, degeneracy_threshold)
        
        num_state_virt = [i + 1 for i in range(len(virtual_values))]

        # Plot occupied states
        for energy, count, states in zip(unique_energies, degenerate_counts, state_numbers):
            if count == 1:
                ax.hlines(energy, xmin=0, xmax=1, color='blue', linewidth=1.5)
                ax.text(0.48, energy + 0.075, f'{states[0]}', va='center', clip_on=True, fontsize=8)
            else:
                ax.hlines(energy, xmin=0, xmax=0.45, color='blue', linewidth=1.5)
                ax.hlines(energy, xmin=0.55, xmax=1, color='blue', linewidth=1.5)
                for j, state in enumerate(states):
                    offset = (0.25 if j % 2 == 0 else 0.75)
                    ax.text(offset, energy + 0.075, f'{state}', va='center', ha='center', clip_on=True, fontsize=8)

        # Plot virtual states
        for energy, num_state in zip(virtual_values, num_state_virt):
            ax.hlines(energy, alpha=1, linewidth=1.5, xmin=0, xmax=1, color='firebrick')
            ax.text(0.48, energy + 0.075, f'{num_state}', va='center', clip_on=True, fontsize=8)

def create_energy_spectra(occupied_values_list, virtual_values_list, titles, split_alpha_beta_list=None, 
                          alpha_occupied_list=None, beta_occupied_list=None, 
                          alpha_virtual_list=None, beta_virtual_list=None, title=None, figsize=None,
                          save_path=None, num_plots=None):

    if figsize and num_plots:
        fig, axes = plt.subplots(num_plots[0], num_plots[1], figsize=(figsize[0] * len(occupied_values_list), figsize[1]), sharey=True)
    else:
        fig, axes = plt.subplots(1, len(occupied_values_list), figsize=(5 * len(occupied_values_list), 6), sharey=True)
    plt.style.use('science')

    for i, ax in enumerate(axes):
        if split_alpha_beta_list and split_alpha_beta_list[i]:
            plot_energy_spectrum(ax, occupied_values_list[i], virtual_values_list[i], titles[i], 
                                 split_alpha_beta=True, alpha_occupied=alpha_occupied_list[i], 
                                 beta_occupied=beta_occupied_list[i], 
                                 alpha_virtual=alpha_virtual_list[i], 
                                 beta_virtual=beta_virtual_list[i])
        else:
            plot_energy_spectrum(ax, occupied_values_list[i], virtual_values_list[i], titles[i])

    fig.text(0, 0.5, r'Orbital energies (eV)', va='center', rotation='vertical', fontsize=15)
    if title:
        fig.text(0.43, 1, title, va='center', fontsize=15)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)

    plt.show()

# Example usage
if __name__ == '__main__':
    # Example data
    occupied_values_NMI2 = [-7.0, -6.8, -6.5, -6.2, -5.9, -5.7]
    virtual_values_NMI2 = [-2.5, -2.0, -1.5]
    occupied_values_singlet_GS = [-7.0, -6.9, -6.7, -6.3, -6.0]
    virtual_values_singlet_GS = [-2.5, -2.0, -1.8]
    alpha_occupied_values_triplet = [-7.1, -6.7, -6.3, -6.0]
    beta_occupied_values_triplet = [-7.0, -6.8, -6.4, -6.1]
    alpha_virtual_values_triplet = [-2.7, -2.2, -1.9]
    beta_virtual_values_triplet = [-2.6, -2.1, -1.8]

    occupied_values_list = [
        occupied_values_NMI2,
        occupied_values_singlet_GS,
        alpha_occupied_values_triplet  # For third plot, will be used in split mode
    ]
    virtual_values_list = [
        virtual_values_NMI2,
        virtual_values_singlet_GS,
        alpha_virtual_values_triplet  # For third plot, will be used in split mode
    ]
    titles = ['NMI2', 'NMI Singlet G-S', 'NMI Triplet']
    split_alpha_beta_list = [False, False, True]  # Enable splitting only for the third plot
    alpha_occupied_list = [None, None, alpha_occupied_values_triplet]
    beta_occupied_list = [None, None, beta_occupied_values_triplet]
    alpha_virtual_list = [None, None, alpha_virtual_values_triplet]
    beta_virtual_list = [None, None, beta_virtual_values_triplet]

    create_energy_spectra(occupied_values_list, virtual_values_list, titles, 
                          split_alpha_beta_list, alpha_occupied_list, beta_occupied_list, 
                          alpha_virtual_list, beta_virtual_list, 
                          save_path='energy_spectra.pdf')

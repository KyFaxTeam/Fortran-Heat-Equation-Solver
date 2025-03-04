#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
import glob
import os
import re
import argparse

def extract_step_number(filename):
    """Extract step number from filename."""
    match = re.search(r'heat_solution_(\d+)\.dat', os.path.basename(filename))
    if match:
        return int(match.group(1))
    return 0

def ensure_directories_exist():
    """Create necessary directories if they don't exist."""
    # Ensure viz directory exists
    viz_dir = 'viz'
    if not os.path.exists(viz_dir):
        os.makedirs(viz_dir)
        print(f"Created directory '{viz_dir}' for visualizations")

def main():
    parser = argparse.ArgumentParser(description='Plot heat equation solutions.')
    parser.add_argument('-a', '--animate', action='store_true', help='Create animation of all solutions')
    parser.add_argument('-s', '--step', type=int, help='Plot specific step')
    parser.add_argument('-m', '--multi', action='store_true', help='Show multiple time steps on same plot')
    parser.add_argument('-n', '--num_steps', type=int, default=5, help='Number of time steps to show (with --multi)')
    parser.add_argument('-f', '--format', default='gif', choices=['gif', 'mp4'], help='Animation format (gif or mp4)')
    parser.add_argument('-o', '--output', default='animation', help='Output filename (without extension)')
    parser.add_argument('-d', '--dpi', type=int, default=150, help='DPI for saved animations')
    parser.add_argument('-c', '--combined', action='store_true', help='Create a combined heatmap view of all steps')
    args = parser.parse_args()

    # Ensure necessary directories exist
    ensure_directories_exist()

    # Ensure results directory exists
    results_dir = 'results'
    if not os.path.exists(results_dir):
        print(f"Directory '{results_dir}' not found! Make sure to run the simulation first.")
        return

    # Get all solution files from results directory
    solution_files = sorted(glob.glob(f'{results_dir}/heat_solution_*.dat'), key=extract_step_number)
    
    if not solution_files:
        print(f"No solution files found in '{results_dir}' directory!")
        return

    # Estimate time step from parameters_mod.f90
    dt = 0.001  # Default value
    try:
        with open('parameters_mod.f90', 'r') as f:
            for line in f:
                if 'dt =' in line:
                    dt_match = re.search(r'dt\s*=\s*([\d.]+)', line)
                    if dt_match:
                        dt = float(dt_match.group(1))
                    break
    except:
        print("Couldn't extract dt from parameters_mod.f90, using default value.")

    if args.step:
        target_file = f'{results_dir}/heat_solution_{args.step:05d}.dat'
        if os.path.exists(target_file):
            plot_solution_file(target_file, dt)
        else:
            print(f"File {target_file} not found!")
        return
    
    if args.animate:
        save_animation(solution_files, dt, args.format, args.output, args.dpi)
    elif args.multi:
        plot_multiple_steps(solution_files, args.num_steps, dt)
    elif args.combined:
        create_combined_heatmap(solution_files, dt)
    else:
        # Plot just the last solution
        plot_solution_file(solution_files[-1], dt)

def plot_solution_file(filename, dt=0.001):
    """Plot a single solution file."""
    data = np.loadtxt(filename, comments='#')
    x = data[:, 0]
    u = data[:, 1]
    
    step = extract_step_number(filename)
    time = step * dt
    
    plt.figure(figsize=(10, 6))
    plt.plot(x, u, 'b-', linewidth=2)
    plt.title(f'Solution de l\'équation de la chaleur à t = {time:.4f}')
    plt.xlabel('Position (x)')
    plt.ylabel('Température (u)')
    plt.grid(True)

    # Save and display
    output_filename = f'viz/plot_step_{step:05d}.png'
    plt.savefig(output_filename, dpi=150)
    print(f"Plot saved to {output_filename}")


def plot_multiple_steps(filenames, num_steps, dt=0.001):
    """Plot multiple time steps on the same graph."""
    plt.figure(figsize=(12, 8))
    
    # Select evenly spaced files if there are more than num_steps files
    if len(filenames) > num_steps:
        indices = np.linspace(0, len(filenames)-1, num_steps, dtype=int)
        selected_files = [filenames[i] for i in indices]
    else:
        selected_files = filenames
    
    # Colors for different curves
    colors = plt.cm.viridis(np.linspace(0, 1, len(selected_files)))
    
    # Store times for colorbar
    times = []
    
    # Plot each solution with a different color
    for i, filename in enumerate(selected_files):
        data = np.loadtxt(filename, comments='#')
        x = data[:, 0]
        u = data[:, 1]
        
        step = extract_step_number(filename)
        time = step * dt
        times.append(time)
        
        # Don't add individual labels - we'll use a colorbar instead
        plt.plot(x, u, '-', color=colors[i], linewidth=2)
    
    plt.title('Évolution de la solution de l\'équation de la chaleur')
    plt.xlabel('Position (x)')
    plt.ylabel('Température (u)')
    plt.grid(True)
    
    # Create a custom colorbar for time instead of individual legend entries
    sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=min(times), vmax=max(times)))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=plt.gca(), label='Temps (t)')
    
    # Add time range info text
    plt.figtext(0.5, 0.01, f"Temps: {min(times):.4f} à {max(times):.4f}", 
               ha="center", fontsize=10)
    
    # Save and display
    output_filename = 'viz/heat_evolution.png'
    plt.savefig(output_filename, dpi=150)
    print(f"Multi-step plot saved to {output_filename}")
  

def save_animation(filenames, dt=0.001, format='gif', output_name='animation', dpi=150):
    """Create and save an animation of all solutions."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # First frame to setup the plot
    data = np.loadtxt(filenames[0], comments='#')
    x = data[:, 0]
    u = data[:, 1]
    
    line, = ax.plot(x, u, 'b-', linewidth=2)
    
    # Find global min and max for consistent y-axis limits
    y_min = float('inf')
    y_max = float('-inf')
    for filename in filenames:
        data = np.loadtxt(filename, comments='#')
        u = data[:, 1]
        y_min = min(y_min, np.min(u))
        y_max = max(y_max, np.max(u))
    
    # Add some padding to the limits
    padding = 0.1 * (y_max - y_min)
    y_min -= padding
    y_max += padding
    
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel('Position (x)')
    ax.set_ylabel('Température (u)')
    ax.grid(True)
    
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
    
    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text
    
    def update(frame):
        data = np.loadtxt(filenames[frame], comments='#')
        x = data[:, 0]
        u = data[:, 1]
        line.set_data(x, u)
        
        step = extract_step_number(filenames[frame])
        time = step * dt
        time_text.set_text(f't = {time:.4f}')
        
        return line, time_text
    
    ani = animation.FuncAnimation(fig, update, frames=len(filenames),
                                 init_func=init, blit=True, interval=100)
    
    if format == 'gif':
        ani.save(f'viz/{output_name}.gif', writer='pillow', dpi=dpi)
        print(f"Animation saved as viz/{output_name}.gif")
    else:
        ani.save(f'viz/{output_name}.mp4', writer='ffmpeg', dpi=dpi)
        print(f"Animation saved as viz/{output_name}.mp4")
    
    plt.close()

def create_combined_heatmap(filenames, dt=0.001):
    """Create a heatmap visualization of the solution over time."""
    # Load all data
    all_u = []
    times = []
    
    for filename in filenames:
        data = np.loadtxt(filename, comments='#')
        all_u.append(data[:, 1])
        
        step = extract_step_number(filename)
        times.append(step * dt)
    
    # Convert to a numpy array for plotting
    all_u = np.array(all_u)
    x = data[:, 0]  # Assuming all files have the same x values
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create a heatmap
    im = ax.imshow(all_u, aspect='auto', origin='lower', 
                  extent=[x.min(), x.max(), times[0], times[-1]],
                  cmap='viridis', interpolation='nearest')
    
    # Add colorbar and labels
    cbar = fig.colorbar(im, ax=ax, label='Température (u)')
    ax.set_xlabel('Position (x)')
    ax.set_ylabel('Temps (t)')
    ax.set_title('Évolution spatio-temporelle de la température')
    
    # Save and show
    plt.savefig('viz/heat_heatmap.png', dpi=150)
    print("Heatmap saved as viz/heat_heatmap.png")


if __name__ == "__main__":
    main()

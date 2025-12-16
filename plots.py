import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import


def load_stats(path):
    data = np.genfromtxt(path, delimiter=",", names=True)
    if data.ndim == 0:
        data = data.reshape(1)
    return data


def plot_time_series(data, output_dir):
    time = data["time"]
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

    axes[0].plot(time, data["pmf_energy"], label="PMF energy")
    axes[0].set_ylabel("PMF energy")
    axes[0].legend()

    axes[1].plot(time, data["pmf_force"], label="PMF force", color="tab:orange")
    axes[1].plot(time, data["center_distance"], label="Chain COM distance", color="tab:green")
    axes[1].set_ylabel("Force / Distance")
    axes[1].legend()

    axes[2].plot(time, data["lj_energy"], label="LJ energy")
    axes[2].plot(time, data["spring_energy"], label="Spring energy")
    axes[2].plot(time, data["rg_a"], label="Rg chain A", linestyle=":")
    axes[2].plot(time, data["rg_b"], label="Rg chain B", linestyle=":")
    axes[2].set_ylabel("Energy / Length")
    axes[2].set_xlabel("Time")
    axes[2].legend()

    fig.tight_layout()
    path = output_dir / "summary_timeseries.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def read_xyz(path):
    frames = []
    with open(path, "r") as fh:
        while True:
            count_line = fh.readline()
            if not count_line:
                break
            atom_count = int(count_line.strip())
            meta = fh.readline().strip()
            positions = []
            for _ in range(atom_count):
                parts = fh.readline().split()
                label = parts[0]
                coords = [float(val) for val in parts[1:4]]
                positions.append((label, coords))
            frames.append((meta, positions))
    return frames


def animate_xyz(frames, output_dir, movie_name="trajectory.gif"):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d")

    def init():
        ax.clear()
        ax.set_xlim3d(-20, 20)
        ax.set_ylim3d(-5, 5)
        ax.set_zlim3d(0, 40)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        return []

    def update(frame):
        ax.clear()
        meta, positions = frame
        xs_a, ys_a, zs_a = [], [], []
        xs_b, ys_b, zs_b = [], [], []
        for label, coords in positions:
            if label == "A":
                xs_a.append(coords[0]); ys_a.append(coords[1]); zs_a.append(coords[2])
            else:
                xs_b.append(coords[0]); ys_b.append(coords[1]); zs_b.append(coords[2])
        ax.scatter(xs_a, ys_a, zs_a, color="tab:blue", s=10, label="Chain A")
        ax.scatter(xs_b, ys_b, zs_b, color="tab:red", s=10, label="Chain B")
        ax.set_title(meta)
        ax.legend(loc="upper right")
        return []

    anim = animation.FuncAnimation(fig, update, frames=frames, init_func=init, interval=100, blit=False)
    output_path = output_dir / movie_name
    anim.save(output_path, writer="pillow", fps=10)
    plt.close(fig)
    return output_path


def main():
    parser = argparse.ArgumentParser(description="Plot Brownian dynamics results")
    parser.add_argument("stats", type=Path, help="Path to stats CSV produced by the simulator")
    parser.add_argument("trajectory", type=Path, help="Path to trajectory XYZ file")
    parser.add_argument("--output-dir", type=Path, default=Path("figures"), help="Directory to store generated plots")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    stats_data = load_stats(args.stats)
    ts_path = plot_time_series(stats_data, args.output_dir)
    print(f"Saved summary time series to {ts_path}")

    if args.trajectory.exists():
        frames = read_xyz(args.trajectory)
        if frames:
            movie_path = animate_xyz(frames, args.output_dir)
            print(f"Saved trajectory animation to {movie_path}")
        else:
            print("No frames found in trajectory; skipping animation")
    else:
        print("Trajectory file missing; only generated plots from stats")


if __name__ == "__main__":
    main()

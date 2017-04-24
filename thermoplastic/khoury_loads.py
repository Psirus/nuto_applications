import multiprocessing
import subprocess


def run_experiment(load):
    run = subprocess.Popen(['./khoury_cooling', '../meshes/3D/CylinderFine.msh', str(load)])
    run.wait()


if __name__ == "__main__":
    loads = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0]
    pool = multiprocessing.Pool(4)
    pool.map(run_experiment, loads)

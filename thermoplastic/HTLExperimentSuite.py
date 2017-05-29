#!/usr/bin/env python3
import multiprocessing
import matplotlib.pyplot as plt
import subprocess
import argparse
import os
import numpy as np

Temperatures = np.linspace(0.0, 1000.0, 11)


def get_strength(directory):
    filename = os.path.join(directory, "TopForce.dat")
    force = np.loadtxt(filename)
    area = np.pi * 31.0**2
    return np.max(np.abs(force)) / area


def plot_strengths(directory):
    strengths = []
    for temperature in Temperatures:
        resultDir = os.path.join(directory, str(temperature))
        strengths.append(get_strength(resultDir))

    plt.plot(Temperatures, strengths)
    figureName = os.path.join(directory, "absolute.pdf")
    plt.savefig(figureName)
    
    plt.cla()
    strengths /= strengths[0]
    plt.plot(Temperatures, strengths)
    figureName = os.path.join(directory, "relative.pdf")
    plt.savefig(figureName)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", help="The program to execute")
    parser.add_argument("mesh", help="The mesh to run it on")
    parser.add_argument("resultDir", help="Where to write the output to")
    args = parser.parse_args()
    print(args.executable)
    print(args.mesh)
    print(args.resultDir)
    os.mkdir(args.resultDir)

    def run_experiment(temperature):
        resultDir = os.path.join(args.resultDir, str(temperature))
        os.mkdir(resultDir)
        run = subprocess.Popen([args.executable, args.mesh, resultDir, str(temperature)])
        print(run.args)
        run.wait()

    pool = multiprocessing.Pool(2)
    pool.map(run_experiment, Temperatures)

import os
import time

import mendeleev as md

import psi4
import numpy as np

import matplotlib.pyplot as plt

frameNumber = 0


def read_atom_list(filename):
    atom_list = []

    with open(filename) as file:
        for line in file:
            x, y, z, atom = line.split(",")
            # *10 factor to convert units
            atom_list.append((float(x) * 10, float(y) * 10, float(z) * 10, md.element(int(atom)).symbol))

    global frameNumber

    # dont show plot for too small atom lists
    if len(atom_list) <= 3:
        return atom_list

    atom_colors = {
        "H": "#0000ff",
        "He": "#ff33ee",
        "N": "#ff00ff",
        "O": "#ff0000",
        "Ne": "#00ff33"
    }

    xs = [x for (x, y, z, a) in atom_list]
    ys = [y for (x, y, z, a) in atom_list]
    zs = [z for (x, y, z, a) in atom_list]

    colors = [atom_colors[a] for (x, y, z, a) in atom_list]

    ax = plt.axes(projection="3d")

    ax.scatter(xs, ys, zs, c=colors)

    plt.savefig(f"../animation/frame{frameNumber}.png")
    frameNumber += 1

    plt.pause(1)

    return atom_list


def save_gradient(filename, gradient, energy):
    n, m = np.shape(gradient)

    g = 49614.7526 * gradient

    with open(filename, "w") as file:
        file.write(f"{energy}\n")

        for i in range(n):
            file.write(f"{g[(i, 0)]}, {g[(i, 1)]}, {g[(i, 2)]}\n")


def save_energy(filename, energy):
    with open(filename, "w") as file:
        file.write(f"{energy * psi4.constants.hartree2kJmol}")


def hartree_fock(atom_list):
    geometry = ""

    for x, y, z, atom in atom_list:
        atom_string = f"{atom} {x} {y} {z}\n"
        geometry += atom_string

    psi4.geometry(geometry)
    psi4.set_options({'basis': '3-21g', "freeze_core": "true"})

    return psi4.gradient('scf'), psi4.variable("scf total energy")


def save_atom_information(filename, atom_types):
    with open(filename, "w") as file:
        for atom in atom_types:
            elem = md.element(atom)
            file.write(f"{elem.atomic_number} {elem.symbol} {elem.atomic_weight} \n")


def optimize_geometry(atom_list):
    geometry = ""

    for x, y, z, atom in atom_list:
        atom_string = f"{atom} {x} {y} {z}\n"
        geometry += atom_string

    psi4.geometry(geometry)
    psi4.set_options({'basis': '3-21g', "freeze_core": "true"})

    energy = psi4.optimize("scf")

    print("Optimized energy is ", energy)

    exit()


def delete_old_files():
    animation_path = "../animation/"
    for file in os.listdir(animation_path):
        if ".png" in file:
            os.remove(animation_path + file)


def main():

    psi4.set_memory("15 GB")

    psi4.set_num_threads(8)
    # psi4.core.be_quiet() #Disable output

    delete_old_files()

    position_file_path = "../tmp/positions.txt"
    gradient_file_path = "../tmp/gradient.txt"
    atom_info_file_path = "../tmp/atom_info.txt"

    while True:
        print("Waiting x file...")

        # Wait until position file exists
        while not os.path.isfile(position_file_path):
            time.sleep(0.1)

        atom_list = read_atom_list(position_file_path)

        os.remove(position_file_path)

        # optimize_geometry(atom_list)  # Optimize Geometry was used to get the relaxed molecules

        atom_types = {atom for (_x, _y, _z, atom) in atom_list}

        print("Calculating...")
        gradient, energy = hartree_fock(atom_list)
        print(energy)

        save_gradient(gradient_file_path, gradient.np, energy)
        save_atom_information(atom_info_file_path, atom_types)


if __name__ == "__main__":
    main()

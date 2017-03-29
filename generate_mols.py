#!/usr/bin/env python2
#
# MIT License
#
# Copyright (c) 2017 Anders Steen Christensen
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import cPickle
import fml

def get_energies(filename):
    """ Returns a dictionary with heats of formation for each xyz-file.
    """

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    energies = dict()

    for line in lines:
        tokens = line.split()

        xyz_name = tokens[0]
        hof = float(tokens[1])

        energies[xyz_name] = hof

    return energies

def get_shieldings(filename):
    """ Returns a dictionary with nmr shieldings for each xyz-file.
    """

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    shieldings = dict()

    for line in lines:
        xyz_name = line[:8]
        shielding = eval(line[8:])

        shieldings[xyz_name] = shielding

    return shieldings

if __name__ == "__main__":

    # Parse file containing PBE0/def2-TZVP heats of formation and xyz filenames
    data = get_energies("hof_qm7.txt")

    # Generate a list of fml.Molecule() objects
    energy_mols = []

    for xyz_file in sorted(data.keys()):

        print xyz_file, data[xyz_file]

        # Initialize the fml.Molecule() objects
        mol = fml.Molecule()
        mol.read_xyz("qm7_xyz/" + xyz_file)

        # Associate a property (heat of formation) with the object
        mol.properties = data[xyz_file]

        energy_mols.append(mol)

    pickle_filename = "mols_hof.cpickle"

    # Save molecules as cPickle
    with open(pickle_filename, "wb") as f:
        cPickle.dump(energy_mols, f, protocol=2)

    print "Wrote", pickle_filename


    # Parse file containing PBE0/def2-TZVP shieldings and xyz filenames
    data = get_shieldings("nmr_qm7.txt")

    # Generate a list of fml.Molecule() objects
    shielding_mols = []

    for xyz_file in sorted(data.keys()):

        print xyz_file, data[xyz_file]

        # Initialize the fml.Molecule() objects
        mol = fml.Molecule()
        mol.read_xyz("qm7_xyz/" + xyz_file)

        # Associate a property (NMR shieldings) with the object
        mol.properties = data[xyz_file]

        shielding_mols.append(mol)

    pickle_filename = "mols_nmr.cpickle"

    # Save molecules as cPickle
    with open(pickle_filename, "wb") as f:
        cPickle.dump(shielding_mols, f, protocol=2)

    print "Wrote", pickle_filename

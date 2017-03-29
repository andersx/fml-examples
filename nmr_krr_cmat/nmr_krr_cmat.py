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

from __future__ import print_function

import sys
import time
import cPickle
import numpy as np
import fml
from fml.kernels import gaussian_kernel
from fml.math import cho_solve, l2_distance

if __name__ == "__main__":

    # Load molecules from cpickle
    pickle_filename = "../mols_nmr.cpickle"
    with open(pickle_filename, "rb") as f:
        mols = cPickle.load(f)

    # Generate descriptor for each molecule
    for mol in mols:

        # This is a list of Coulomb matrices sorted by distance
        # to each query atom.
        mol.generate_atomic_coulomb_matrix()

    # Shuffle molecules
    np.random.seed(666)
    np.random.shuffle(mols)

    # Make training and test sets
    n_test  = 1000
    n_train = 4000
    
    Xall = []
    Yall = []

    target_type = "H"

    for mol in mols:
        for i, atomtype in enumerate(mol.atomtypes):
            if atomtype == target_type:
                Xall.append(mol.atomic_coulomb_matrix[i])
                Yall.append(mol.properties[i])

    # Vectors of descriptors for training and test sets - note transposed
    # for enhanced speed in kernel evaluation
    X = np.array(Xall[:n_train]).T
    Xs = np.array(Xall[-n_test:]).T

    # Vectors of properties for training and test sets
    Y = np.array(Yall[:n_train])
    Ys = np.array(Yall[-n_test:])

    # Set hyper-parameters
    sigma = 10**(4.2)
    llambda = 10**(-10.0)

    # print "Calculating K-matrix           ...",
    print(u"Calculating: K\u1D62\u2C7C = k(q\u1D62, q\u2C7C)          ... ", end="")
    sys.stdout.flush() 
    start = time.time()

    # Gaussian kernel usually better for atomic properties
    # K = gaussian_kernel(X, X, sigma)

    # Alternatively, just calculate the L2 distance, and convert
    # to kernel matrix (e.g. for multiple sigmas, etc):
    D = l2_distance(X,X)
    D /= -2.0 * sigma * sigma
    D = np.exp(K)

    print ("%7.2f seconds" % (time.time() - start) )

    for i in xrange(n_train):
        K[i,i] += llambda

    print( u"Calculating: \u03B1 = (K + \u03BBI)\u207B\u00B9 y         ... ", end="")
    sys.stdout.flush()
    start = time.time()
    alpha = cho_solve(K,Y)
    print ("%7.2f seconds" % (time.time() - start) )

    print(u"Calculating: K*\u1D62\u2C7C = k(q\u1D62, q*\u2C7C)        ... ", end="")
    sys.stdout.flush()
    start = time.time()
    Ks = gaussian_kernel(X, Xs, sigma)
    print ("%7.2f seconds" % (time.time() - start) )


    print( u"Calculating: y* = (K*)\u1D40\u00B7\u03B1             ... ", end="")
    sys.stdout.flush()
    start = time.time()
    Y_tilde = np.dot(Ks.transpose(), alpha)
    print ("%7.2f seconds" % (time.time() - start) )
    
    rmsd = np.sqrt(np.mean(np.square(Ys - Y_tilde)))
    print("RMSD = %6.2f ppm" % rmsd)

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
from fml.kernels import laplacian_kernel
from fml.math import cho_solve

if __name__ == "__main__":

    # Load molecules from cpickle
    pickle_filename = "../mols_hof.cpickle"
    with open(pickle_filename, "rb") as f:
        mols = cPickle.load(f)

    # Generate descriptor for each molecule
    for mol in mols:
        mol.generate_coulomb_matrix()

    # Shuffle molecules
    np.random.seed(666)
    np.random.shuffle(mols)

    # Make training and test sets
    n_test  = 1000
    n_train = 1000
    
    training = mols[:n_train]
    test  = mols[-n_test:]

    # List of descriptors for training set -- note transposed, because quirk
    X  = np.array([mol.coulomb_matrix for mol in training]).T

    # List of properties for training set
    Y = np.array([mol.properties for mol in training])

    # List of descriptors for test set -- note transposed, because quirk
    Xs = np.array([mol.coulomb_matrix for mol in test]).T

    # List of properties for test set
    Ys = np.array([mol.properties for mol in test])

    # Set hyper-parameters
    sigma = 10**(4.2)
    llambda = 10**(-10.0)

    # print "Calculating K-matrix           ...",
    print(u"Calculating: K\u1D62\u2C7C = k(q\u1D62, q\u2C7C)          ... ", end="")
    sys.stdout.flush() 
    start = time.time()
    K = laplacian_kernel(X, X, sigma)
    print ("%7.2f seconds" % (time.time() - start) )

    for i in xrange(len(training)):
        K[i,i] += llambda

    print( u"Calculating: \u03B1 = (K + \u03BBI)\u207B\u00B9 y         ... ", end="")
    sys.stdout.flush()
    start = time.time()
    alpha = cho_solve(K,Y)
    print ("%7.2f seconds" % (time.time() - start) )

    print(u"Calculating: K*\u1D62\u2C7C = k(q\u1D62, q*\u2C7C)        ... ", end="")
    sys.stdout.flush()
    start = time.time()
    Ks = laplacian_kernel(X, Xs, sigma)
    print ("%7.2f seconds" % (time.time() - start) )


    print( u"Calculating: y* = (K*)\u1D40\u00B7\u03B1             ... ", end="")
    sys.stdout.flush()
    start = time.time()
    Y_tilde = np.dot(Ks.transpose(), alpha)
    print ("%7.2f seconds" % (time.time() - start) )
    
    rmsd = np.sqrt(np.mean(np.square(Ys - Y_tilde)))
    print("RMSD = %6.2f kcal/mol" % rmsd)

Probabilistic Nordsieck method
==============================

About
-----

This software release supplements the [paper](https://arxiv.org/abs/1610.05261):

    Michael Schober, Simo Särkkä, Philipp Hennig:
    "A probabilistic model for the numerical solution of initial value problems", 
    2017.

The numerical implementation is provided in the Matlab programming
environment. The probabilistic Nordsieck method is implemented in the
function

    odeFilter

whose interface resembles other numerical differential equation
solvers available in Matlab. Other functions in

    solver/filter/

implement additional functionality associated with the filter output,
such as computing the smoothing distribution and sampling from the
predictive posterior.

To reproduce the illustrative plots from Section 2 of the paper, open
a Matlab instance and type

    setup; Sec2Figure

Code to reproduce the benchmark comparison is partially copyrighted by
Matlab and, thus, cannot be publicly released. We are working to
provide an alternative soon.

Further symbolic algebra code written in Python is provided to
check the derivations.

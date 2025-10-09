A Binder repo for the stuff I made with SageMath
================================================

If you click the button below, an environment in Binder will be generated to interact with my Jupyter Notebooks.

.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/joelcrey/JupyterNotebooks/HEAD?urlpath=tree


My Jupyter Notebooks
____________________
1. The notebook 'hanmonskyalgorithm.ipynb' contains just one cell. The main functions are HK_diagonal and eHK_diagonal which compute the Hilbert–Kunz function and Hilbert–Kunz multiplicity of a diagonal hypersurface via the Han–Monsky algorithm. You just need to introduce the characteristic and the collection of exponents. For example: to compute the Hilbert–Kunz function of x^3+y^3+z^3 in characteristic 5, simply introduce HK_diagonal(5,[3,3,3]), which will output 'HK_(e)=9/4*5^(2*e)-5/4 for e>=0'.
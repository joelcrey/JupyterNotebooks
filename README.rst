A Binder repo for the stuff I made with SageMath
================================================

These are Jupyter notebooks of Sagemath code to compute Hilbert–Kunz functions, multiplicities and F-signatures of diagonal hypersurfaces via Han–Monsky and Chiang–Hung algorithm.

If you cannot run it locally, click the button below to open a Binder environment to interact with these Jupyter Notebooks.

.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/joelcrey/JupyterNotebooks/HEAD?urlpath=tree


My Jupyter Notebooks
____________________
1. The notebook 'Han_Monsky_Algorithm.ipynb' allows you to compute the Hilbert–Kunz function, Hilbert–Kunz multiplicity and F-signature of any hypersurface of the form x_1^d_1+...+x_n^d_n for any fixed characteristic p. The main functions are full_HK_function_diagonal(p,[d_1,...,d_n]), eHK_diagonal(p,[d_1,...,d_n]) and s_diagonal(p,[d_1,...,d_n]), which compute the Hilbert–Kunz function, Hilbert–Kunz multiplicity and F-signature of said diagonal hypersurfaces via the Han–Monsky algorithm [HM93]. I originally implemented the algorithm from Han and Monsky's paper [HM93] and later asked Claude to improve the notebook's presentation; part 1b is entirely Claude's. The code has some worked-out examples at the end.

2. The notebook 'Chiang_Hung_Algorithm_by_Claude.ipynb' (according to Claude) allows you to compute the Hilbert–Kunz function and multiplicity of any hypersurface of the form (x_(1,1)...x_(1,d_1))^e_1+...+(x_(m,1)...x_(m,d_m))^e_m (i.e., sum of monomials of disjoint variables); in particular, diagonal hypersurfaces are an example of this, but this algorithm is slower than Han and Monsky's one. Disclaimer: The function in this notebook was not implemented by me; it was implemented entirely by Claude after being prompted to implement the algorithm developed by Chiang and Hung in their 1996 paper [CH96], so I claim no correctness of the algorithm; however, it checks out in all cases from [CH96] I could try, and also with numerical experiments for first values of the Hilbert–Kunz functions of several other hypersurfaces.

References
____________________
[HM93] Han, C., & Monsky, P. (1993). Some surprising Hilbert-Kunz functions. Mathematische Zeitschrift, 214(1), 119-135.
[CH96] Chiang, L., & Hung, Y. C. (1998). On Hilbert-Kunz functions of some hypersurfaces. Journal of Algebra, 199(2), 499-527.


A Binder repo for the stuff I made with SageMath
================================================

These are Jupyter notebooks of Sagemath code to compute Hilbert–Kunz functions, multiplicities and F-signatures of diagonal hypersurfaces via Han–Monsky and Chiang–Hung algorithm.

If you cannot run it locally, click the button below to open a Binder environment to interact with these Jupyter Notebooks.

.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/joelcrey/JupyterNotebooks/HEAD?urlpath=tree


My Jupyter Notebooks
____________________
All the notebooks are in the folder 'notebooks'. The notebooks 01, 02 and 03 use the module 'hilbert_kunz.py', which must be in the same folder. They follow the notes 'han_monsky.pdf' [Notes].

1. The notebook 'hanmonskyalgorithm.ipynb' is my original implementation of the algorithm from Han and Monsky's paper [HM93], before asking Claude to improve it. It computes the Hilbert–Kunz multiplicity of any diagonal hypersurface (eHK_diagonal(p,[e_0,...,e_d])) and the Hilbert–Kunz function of any diagonal hypersurface (HK_diagonal(p,[e_0,...,e_d])). It usually does not work if the characteristic or the exponents of the diagonal hypersurface are large.

2. The notebook '01_representation_ring.ipynb' allows you to compute in the Han–Monsky ring (or representation ring) [HM93]: products of deltas and lambdas, changes of basis between them, and Yoshida's matrix representation.

3. The notebook '02_han_monsky.ipynb' allows you to compute the Hilbert–Kunz function, Hilbert–Kunz multiplicity, F-signature function and F-signature of any hypersurface of the form x_1^d_1+...+x_n^d_n for any fixed characteristic p, via the Han–Monsky algorithm [HM93]. The main functions are HK_function_diagonal(f) and FS_function_diagonal(f), where f is the polynomial, and eHK_diagonal(p,[d_1,...,d_n]) and F_signature_diagonal(p,[d_1,...,d_n]), which compute the Hilbert–Kunz multiplicity and F-signature of said diagonal hypersurfaces. This is a newer version of the notebook above that I made with Claude. The code has some worked-out examples at the beginning, and the algorithm is explained afterwards.

4. The notebook '03_chiang_hung.ipynb' (according to Claude) allows you to compute the Hilbert–Kunz function and multiplicity of any hypersurface of the form (x_(1,1)...x_(1,d_1))^e_1+...+(x_(m,1)...x_(m,d_m))^e_m (i.e., sum of monomials of disjoint variables); in particular, diagonal hypersurfaces are an example of this, but this algorithm is slower than Han and Monsky's one. The main function is HK_function_disjoint_monomials(f). Disclaimer: The function in this notebook was not implemented by me; it was implemented entirely by Claude after being prompted to implement the algorithm developed by Chiang and Hung in their 1998 paper [CH98], so I claim no correctness of the algorithm; however, it checks out in all cases from [CH98] I could try, and also with numerical experiments for first values of the Hilbert–Kunz functions of several other hypersurfaces.

References
____________________
[Notes] Castillo-Rey, J. The Hilbert–Kunz function and F-signature function of the diagonal hypersurfaces. Notes ('notebooks/han_monsky.pdf').
[HM93] Han, C., & Monsky, P. (1993). Some surprising Hilbert-Kunz functions. Mathematische Zeitschrift, 214(1), 119-135.
[CH98] Chiang, L., & Hung, Y. C. (1998). On Hilbert-Kunz functions of some hypersurfaces. Journal of Algebra, 199(2), 499-527.

The full list of references is at the end of each notebook.

# A tenth-degree polynomial equation — revised Beamer slides

Editable 15-slide presentation, reconstructed from the supplied PDF and revised in accordance with the author's mathematical typography and pedagogical sequence. The source uses Beamer in 16:9, Libertinus Serif/Math, and TikZ.

## Revision details

- Standardized the manuscript to American English: **factorization**, **organized**, **recognized**, **generalization**.
- Eq. (2) now stays on one line. Eq. (9) ends with the explicit product `(1+x+x^2)(1+x^4+x^8)`.
- The two geometric-series identities are numbered (10a) and (10b); Eq. (11) retains its numbering.
- Slide 7 now has a larger, non-overlapping root diagram, justified explanatory text, and Eq. (12) left-aligned using `flalign` and an `aligned` system with a left brace.
- The original slide 11 is divided into three pedagogical slides: real factorization, rational factorization, and the definition of monic/cyclotomic polynomials.
- Reference hanging indents were regularized. The final reference is now T. Needham, *Visual Complex Analysis* (Oxford University Press, 1997). All display equations are numbered without mathematical boxes.
- Mathematics was checked symbolically with SymPy; the PDF was compiled twice with XeLaTeX.

The associated earlier presentation can be found at [video 041 in youtube_dr_olsp](https://github.com/ozsp12/youtube_dr_olsp/tree/main/videos/041_math_algebra_polynomial_equation_10th_degree).

## Reproducible build

Required: XeLaTeX, Beamer, fontspec, unicode-math, Libertinus, TikZ, Python and SymPy.

```sh
make
make verify
```

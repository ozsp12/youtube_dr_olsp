.PHONY: all verify clean
all:
	xelatex -interaction=nonstopmode -halt-on-error polynomial_10th_degree.tex
	xelatex -interaction=nonstopmode -halt-on-error polynomial_10th_degree.tex
verify:
	python checks/verify_polynomial.py
clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.synctex.gz

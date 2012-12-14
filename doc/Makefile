RST=$(wildcard *.rst)
PDF=$(patsubst %.rst,%.pdf,$(RST))

all: $(PDF)

clean:
	rm -f *.pdf

.PHONY: clean

%.pdf: %.rst slides.style
	rst2pdf $< -b1 -o $@ -s slides.style


RST=$(wildcard *.rst)
PDF=$(patsubst %.rst,%.docx,$(RST))

all: $(PDF)

clean:
	rm -f *.pdf

.PHONY: clean

%.docx: %.rst slides.style
#	rst2pdf --use-floating-images $< -b1 -o $@ -s slides.style
	pandoc -f rst -t docx $< -o $@

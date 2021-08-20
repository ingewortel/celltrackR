
all : vignette-pages

vignettes=reading-converting-data QC data-QC ana-methods clustering simulation


vignette-pages : $(foreach v, $(vignettes), vignettes-out/$(v).html )


vignettes-out : 
	mkdir -p $@

vignettes-out/%.html : page-header.html navbar.html vignettes/%.html | vignettes-out
	cat $< > vignettes-out/$*-tmp.html && \
	tail -n +12 vignettes/$*.html >> vignettes-out/$*-tmp.html && \
	sed '/<body>/ r navbar.html' vignettes-out/$*-tmp.html > $@

.PHONY: clean cleaner

all:
	python3 support/tutorial_tools/notebook/process_notebook.py index cross-link_ms-ambiguity cross-link_ms

clean:
	rm -f index.md cross-link_ms.md cross-link_ms-ambiguity.md Doxyfile
	rm -f excluded.None.xl.db included.None.xl.db missing.None.xl.db xlinks.csv
	rm -rf html

cleaner:
	rm -f index.ipynb cross-link_ms.ipynb cross-link_ms-ambiguity.ipynb cross-link_ms.py cross-link_ms-ambiguity.py cross-link_ms.md cross-link_ms-ambiguity.md index.md Doxyfile
	rm -f excluded.None.xl.db included.None.xl.db missing.None.xl.db xlinks.csv
	rm -rf html


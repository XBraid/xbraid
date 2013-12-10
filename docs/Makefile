# Requires doxygen version at least 1.8
DOX=doxygen
DOXVERSION := $(shell expr `doxygen --version | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 10800)

clean_latex:
	# sed command cleans up latex file some
	cd user_manual;      sed '/.*8md.*/d' refman.tex > refman.temp.tex; mv refman.temp.tex refman.tex
	cd developer_manual; sed '/.*8md.*/d' refman.tex > refman.temp.tex; mv refman.temp.tex refman.tex

user_manual: user_manual.conf 
ifeq "$(DOXVERSION)" "1"
	$(DOX) user_manual.conf
	make clean_latex
	cd user_manual; make; cp refman.pdf ../user_manual.pdf
else
	@echo "Doxygen version too old"
endif	

developer_manual: developer_manual.conf 
ifeq "$(DOXVERSION)" "1"
	$(DOX) developer_manual.conf
	make clean_latex 
	cd developer_manual; make; cp refman.pdf ../developer_manual.pdf
else
	@echo "Doxygen version too old"
endif	

clean:
	rm -rf developer_manual
	rm -rf user_manual
	rm -f user_manual.pdf
	rm -f developer_manual.pdf


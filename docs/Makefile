#BHEADER**********************************************************************
# Copyright (c) 2013,  Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# This file is part of XBraid.  See file COPYRIGHT for details.
#
# XBraid is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.
#
#EHEADER**********************************************************************

# Requires doxygen version at least 1.8
DOX=doxygen
DOXVERSION := $(shell expr `doxygen --version | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 10800)

clean_latex:
	# sed command cleans up latex file some
	if test -d user_manual; then cd user_manual; sed '/.*8md.*/d' refman.tex > refman.temp.tex; mv refman.temp.tex refman.tex; fi
	if test -d developer_manual; then cd developer_manual; sed '/.*8md.*/d' refman.tex > refman.temp.tex; mv refman.temp.tex refman.tex; fi
	# sed command changes pdflatex to nonstop mode (important for regression testing)
	if test -d user_manual; then cd user_manual; sed 's/pdflatex/pdflatex --interaction=nonstopmode/g' Makefile > Makefile.temp; mv Makefile.temp Makefile; fi
	if test -d developer_manual; then cd developer_manual; sed 's/pdflatex/pdflatex --interaction=nonstopmode/g' Makefile > Makefile.temp; mv Makefile.temp Makefile; fi
	# sed command changes makeindex to quiet mode (important for regression testing)
	if test -d user_manual; then cd user_manual; sed 's/makeindex/makeindex -q/g' Makefile > Makefile.temp; mv Makefile.temp Makefile; fi
	if test -d developer_manual; then cd developer_manual; sed 's/makeindex/makeindex -q/g' Makefile > Makefile.temp; mv Makefile.temp Makefile; fi
	
user_manual: user_manual.conf 
ifeq "$(DOXVERSION)" "1"
	$(DOX) user_manual.conf
	make clean_latex
	cd user_manual; make; cp refman.pdf ../user_manual.pdf
else
	@echo "Error: Doxygen version too old"
endif	

developer_manual: developer_manual.conf 
ifeq "$(DOXVERSION)" "1"
	$(DOX) developer_manual.conf
	make clean_latex 
	cd developer_manual; make; cp refman.pdf ../developer_manual.pdf
else
	@echo "Error: Doxygen version too old"
endif	

clean:
	rm -rf developer_manual
	rm -rf user_manual
	rm -f user_manual.pdf
	rm -f developer_manual.pdf


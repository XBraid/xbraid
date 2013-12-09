# Requires doxygen version at least 1.8
DOX=doxygen
DOXVERSION := $(shell expr `doxygen --version | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 10800)

user_manual: user_manual.conf postprocess_latex
ifeq "$(DOXVERSION)" "1"
	$(DOX) user_manual.conf
	# sed command cleans up latex file some
	sed -i '/.*8md.*/d' user_manual/refman.tex 2> /dev/null	
	cd user_manual; make; cp refman.pdf ../user_manual.pdf
else
	@echo "Doxygen version too old"
endif	

developer_manual: developer_manual.conf postprocess_latex
ifeq "$(DOXVERSION)" "1"
	$(DOX) developer_manual.conf
	# sed command cleans up latex file some
	sed -i '/.*8md.*/d' developer_manual/refman.tex 2> /dev/null
	cd developer_manual; make; cp refman.pdf ../developer_manual.pdf
else
	@echo "Doxygen version too old"
endif	

clean:
	rm -rf developer_manual
	rm -rf user_manual
	rm -f user_manual.pdf
	rm -f developer_manual.pdf


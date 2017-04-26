
export MKFILE_DIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SRC=src
TARGET_BIN=bin
TOOL_PATHS=$(TARGET_BIN)/tools.py
PROGRAM=$(SVE_DIR)/scripts/sve

SVE_DIR=$(shell pwd)
R_PACKAGE=$(SVE_DIR)/$(SRC)/R-package
R_INSTALL_DIR=$(R_PACKAGE)/packages
R_PACKAGE_DEPEN = bzip2-1.0.6 curl-7.47.1 pcre-8.40 xz-5.2.2 zlib-1.2.9

PERL_LIB=$(SVE_DIR)/$(SRC)/perl-lib
PERL_LIB_DEPEN = GD-2.52 GDTextUtil-0.86 GDGraph-histogram-1.1 GDGraph-1.54

TARBALLS = jre1.8.0_51 picard-tools-2.5.0 svtoolkit_2.00.1736 CNVnator_v0.3.3 samtools-0.1.19 breakdancer-1.4.5
SUBDIRS = bwa speedseq htslib samtools samtools-0.1.19 bcftools bedtools2 delly lumpy-sv hydra tigra

# all
all: unzip_tarballs configure build CNVnator_v0.3.3 perl-lib breakdancer
	@test -d $(SVE_DIR)/data || tar -zxvf data.tar.gz # unzip data
	@test -d $(SVE_DIR)/$(TARGET_BIN) || mkdir $(SVE_DIR)/$(TARGET_BIN)
	$(MAKE) tool_paths
	@cp $(PROGRAM) $(SVE_DIR)/$(TARGET_BIN)
.PHONY: all

unzip_tarballs:
	@cd $(SVE_DIR)/$(SRC); \
	for module in $(TARBALLS); do \
		if [ ! -d $$module ]; then \
			echo "- Unzip $$module"; \
			tar -zxvf $$module.tar.gz; \
		fi; \
	done; \
	cd $(SVE_DIR)

configure:
	@echo "- Configuring in htslib"
	@cd $(SVE_DIR)/$(SRC)/htslib && autoheader && autoconf && ./configure --disable-lzma && $(MAKE) 
	@cd $(SVE_DIR)
	@echo "- Configuring in lumpy"
	@cd $(SVE_DIR)/$(SRC)/lumpy-sv/lib/htslib && autoheader && autoconf
	@cd $(SVE_DIR)
	@echo "- Configuring in breakseq2"
	@cd $(SVE_DIR)/$(SRC)/breakseq2 && python setup.py && cd $(SVE_DIR) || true
	@echo "- Configuring in tigra"
	@sed -i "/SAMTOOLS=/d" $(SVE_DIR)/$(SRC)/tigra/Makefile
	@sed -i "/HTSLIB=/d" $(SVE_DIR)/$(SRC)/tigra/Makefile
	@sed -i "/SAMHTSLIB=/d" $(SVE_DIR)/$(SRC)/tigra/Makefile
	@sed -i "2 a SAMHTSLIB=$(SVE_DIR)/$(SRC)/htslib/htslib" $(SVE_DIR)/$(SRC)/tigra/Makefile
	@sed -i "2 a HTSLIB=$(SVE_DIR)/$(SRC)/htslib" $(SVE_DIR)/$(SRC)/tigra/Makefile
	@sed -i "2 a SAMTOOLS=$(SVE_DIR)/$(SRC)/samtools" $(SVE_DIR)/$(SRC)/tigra/Makefile

build:
	@for dir in $(SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $(SRC)/$$dir;\
	done

CNVnator_v0.3.3:
	$(MAKE) --no-print-directory -C $(SRC)/CNVnator_v0.3.3/src/samtools
	$(MAKE) --no-print-directory -C $(SRC)/CNVnator_v0.3.3/src

breakdancer:
	@cd $(SVE_DIR)/$(SRC)/breakdancer-1.4.5; \
	mkdir build; \
	cd build; \
	cmake .. -DCMAKE_BUILD_TYPE=release; \
	$(MAKE)

perl-lib:
	@cd $(PERL_LIB); \
	for module in $(PERL_LIB_DEPEN); do \
		if [ ! -d $$module ]; then \
			echo "- Unzip $$module"; \
			tar -zxvf $$module.tar.gz; \
			cd $(PERL_LIB)/$$module; \
			perl Makefile.PL INSTALL_BASE=$(PERL_LIB); \
			$(MAKE) && $(MAKE) install; \
		fi; \
	done; \

R-package:
	@test -d $(R_INSTALL_DIR) || mkdir -p $(R_INSTALL_DIR)
	@cd $(R_PACKAGE); \
	for module in $(R_PACKAGE_DEPEN); do \
		if [ ! -d $$module ]; then \
			echo "- Unzip $$module"; \
			tar -zxvf $$module.tar.gz; \
		fi; \
	done; \
	test -d R-3.3.3 || tar -zxvf R-3.3.3.tar.gz; \
	cd $(SVE_DIR)
	@sed -i "/CFLAGS=/d" $(R_PACKAGE)/bzip2-1.0.6/Makefile
	@sed -i '23 a CFLAGS=-Wall -Winline -O2 -g -fPIC $$(BIGFILES)' $(R_PACKAGE)/bzip2-1.0.6/Makefile
	@cd $(R_PACKAGE)/bzip2-1.0.6 && $(MAKE) clean && $(MAKE) -f Makefile-libbz2_so && $(MAKE) clean && $(MAKE) && $(MAKE) -n install PREFIX=$(R_INSTALL_DIR) && $(MAKE) install PREFIX=$(R_INSTALL_DIR)
	@cd $(R_PACKAGE)/curl-7.47.1 && ./configure --prefix=$(R_INSTALL_DIR) && $(MAKE) -j3 && $(MAKE) install
	@cd $(R_PACKAGE)/pcre-8.40 && ./configure --prefix=$(R_INSTALL_DIR) --enable-utf8 && $(MAKE) && $(MAKE) install
	@cd $(R_PACKAGE)/xz-5.2.2 && ./configure --prefix=$(R_INSTALL_DIR) && $(MAKE) -j3 && $(MAKE) install
	@cd $(R_PACKAGE)/zlib-1.2.9 && ./configure --prefix=$(R_INSTALL_DIR) && $(MAKE) && $(MAKE) install
	@cd $(R_PACKAGE)/R-3.3.3 && ./configure --prefix=$(R_INSTALL_DIR) LDFLAGS='-L$(R_INSTALL_DIR)/lib' CFLAGS='-I$(R_INSTALL_DIR)/include' && $(MAKE) || true
	@$(MAKE) -C $(R_PACKAGE)/R-3.3.3
	@cd $(SVE_DIR)

tool_paths:
	@echo "TOOLS={}" > $(TOOL_PATHS)
	@echo "TOOLS ['SVE_HOME']      = '$(SVE_DIR)'" >> $(TOOL_PATHS)
	@echo "TOOLS ['BWA']           = '$(SVE_DIR)/$(SRC)/bwa/bwa'" >> $(TOOL_PATHS)
	@echo "TOOLS ['BWA-POSTALT']   = '$(SVE_DIR)/$(SRC)/k8 $(SVE_DIR)/$(SRC)/bwa/bwakit/bwa-postalt.js'" >> $(TOOL_PATHS)
	@echo "TOOLS ['SPEEDSEQ']      = '$(SVE_DIR)/$(SRC)/speedseq/bin/speedseq'" >> $(TOOL_PATHS)
	@echo "TOOLS ['SAMBAMBA']      = '$(SVE_DIR)/$(SRC)/speedseq/bin/sambamba'" >> $(TOOL_PATHS)
	@echo "TOOLS ['SAMTOOLS']      = '$(SVE_DIR)/$(SRC)/samtools/samtools'" >> $(TOOL_PATHS)
	@echo "TOOLS ['TABIX']         = '$(SVE_DIR)/$(SRC)/htslib/tabix'" >> $(TOOL_PATHS)
	@echo "TOOLS ['JAVA-1.8']      = '$(SVE_DIR)/$(SRC)/jre1.8.0_51/bin/java'" >> $(TOOL_PATHS)
	@echo "TOOLS ['PICARD']        = '$(SVE_DIR)/$(SRC)/picard-tools-2.5.0/picard.jar'" >> $(TOOL_PATHS)
	@echo "TOOLS ['CNVNATOR']      = '$(SVE_DIR)/$(SRC)/CNVnator_v0.3.3/src/cnvnator'" >> $(TOOL_PATHS)
	@echo "TOOLS ['CNVNATOR2VCF']  = '$(SVE_DIR)/$(SRC)/CNVnator_v0.3.3/cnvnator2VCF.pl'" >> $(TOOL_PATHS)
	@echo "TOOLS ['DELLY']         = '$(SVE_DIR)/$(SRC)/delly/src/delly'" >> $(TOOL_PATHS)
	@echo "TOOLS ['BCFTOOLS']      = '$(SVE_DIR)/$(SRC)/bcftools/bcftools'" >> $(TOOL_PATHS)
	@echo "TOOLS ['BEDTOOLS']      = '$(SVE_DIR)/$(SRC)/bedtools2/bin/bcftools'" >> $(TOOL_PATHS)
	@echo "TOOLS ['LUMPY-EXPRESS'] = '$(SVE_DIR)/$(SRC)/lumpy-sv/bin/lumpyexpress'" >> $(TOOL_PATHS)
	@echo "TOOLS ['BREAKSEQ']      = '$(SVE_DIR)/$(SRC)/breakseq2/scripts/run_breakseq2.py'" >> $(TOOL_PATHS)
	@echo "TOOLS ['TIGRA']         = '$(SVE_DIR)/$(SRC)/tigra/tigra-sv'" >> $(TOOL_PATHS)
	@echo "TOOLS ['TIGRA-EXT']     = '$(SVE_DIR)/$(SRC)/tigra-ext/TIGRA-ext.pl'" >> $(TOOL_PATHS)
	@echo "TOOLS ['FATO2BIT']      = '$(SVE_DIR)/$(SRC)/faToTwoBit/faToTwoBit_linux'" >> $(TOOL_PATHS)
	@echo "TOOLS ['R_PATH']               = '$(R_PACKAGE)/R-3.3.3'" >> $(TOOL_PATHS)
	@echo "TOOLS ['R_LIB_PATH']               = '$(R_INSTALL_DIR)'" >> $(TOOL_PATHS)
	@echo "TOOLS ['PERL_LIB_PATH']            = '$(PERL_LIB)'" >> $(TOOL_PATHS)
	@echo "TOOLS ['HYDRA_PATH']           = '$(SVE_DIR)/$(SRC)/hydra'" >> $(TOOL_PATHS)
	@echo "TOOLS ['BREAKDANCER_PATH']     = '$(SVE_DIR)/$(SRC)/breakdancer-1.4.5'" >> $(TOOL_PATHS)
	@echo "TOOLS ['TIGRA_PATH']           = '$(SVE_DIR)/$(SRC)/tigra'" >> $(TOOL_PATHS)
	@echo "TOOLS ['BWA_PATH']             = '$(SVE_DIR)/$(SRC)/bwa'" >> $(TOOL_PATHS)
	@echo "TOOLS ['SAMTOOLS_PATH']        = '$(SVE_DIR)/$(SRC)/samtools'" >> $(TOOL_PATHS)
	@echo "TOOLS ['SAMTOOLS-0.1.19_PATH'] = '$(SVE_DIR)/$(SRC)/samtools-0.1.19'" >> $(TOOL_PATHS)
	@echo "TOOLS ['BCFTOOLS_PATH']        = '$(SVE_DIR)/$(SRC)/bcftools'" >> $(TOOL_PATHS)
	@echo "TOOLS ['HTSLIB_PATH']          = '$(SVE_DIR)/$(SRC)/htslib'" >> $(TOOL_PATHS)
	@echo "TOOLS ['JAVA-1.8_PATH']        = '$(SVE_DIR)/$(SRC)/jre1.8.0_51/bin'" >> $(TOOL_PATHS)
	@echo "TOOLS ['GENOME_STRIP_PATH']    = '$(SVE_DIR)/$(SRC)/svtoolkit'" >> $(TOOL_PATHS)
	@echo "TOOLS ['BREAKSEQ_PATH']        = '$(SVE_DIR)/$(SRC)/breakseq2'" >> $(TOOL_PATHS)
	@echo "" >> $(TOOL_PATHS)
	@echo "FILES={}" >> $(TOOL_PATHS)
	@echo "FILES ['GRCH38-EXTRA'] = '$(SVE_DIR)/data/bwakit-GRCh38/hs38DH-extra.fa'"  >> $(TOOL_PATHS)
	@echo "FILES ['GRCH38-ALT'] = '$(SVE_DIR)/data/bwakit-GRCh38/hs38DH.fa.alt'"  >> $(TOOL_PATHS)
	@echo "FILES ['DELLY-HG19'] = '$(SVE_DIR)/$(SRC)/delly/excludeTemplates/human.hg19.excl.tsv'"  >> $(TOOL_PATHS)
	@echo "FILES ['DELLY-HG38'] = '$(SVE_DIR)/$(SRC)/delly/excludeTemplates/human.hg38.excl.tsv'"  >> $(TOOL_PATHS)
	@echo "FILES ['BREAKSEQ-HG19'] = '$(SVE_DIR)/data/breakseq_bplib/breakseq2_bplib_20150129.gff'"  >> $(TOOL_PATHS)
	@echo "FILES ['BREAKSEQ-HG38'] = '$(SVE_DIR)/data/breakseq_bplib/breakseq2_bplib_20150129_hg38.gff'"  >> $(TOOL_PATHS)
clean:
	@echo "- Clean up"
	@for dir in $(SUBDIRS); do \
		if [ -d $$dir ]; then \
			$(MAKE) -C $(SRC)/$$dir clean; \
		fi; \
	done
	@for module in $(TARBALLS); do \
		if [ -d $(SRC)/$$module ]; then \
			rm -rf $(SRC)/$$module; \
		fi; \
	done
	@for module in $(R_PACKAGE_DEPEN); do \
		if [ -d $(R_PACKAGE)/$$module ]; then \
			rm -rf $(R_PACKAGE)/$$module; \
		fi; \
	done;
	@rm -rf $(R_INSTALL_DIR)
	@rm -rf $(R_PACKAGE)/R-3.3.3
	@for module in $(PERL_LIB_DEPEN); do \
		if [ -d $(PERL_LIB)/$$module ]; then \
			rm -rf $(PERL_LIB)/$$module; \
		fi; \
	done;
	@if [ -d $(SVE_DIR)/data ]; then \
		chmod -R 777 $(SVE_DIR)/data; \
		rm -rf $(SVE_DIR)/data; \
	fi
	@rm -rf $(TARGET_BIN)

.PHONY: clean

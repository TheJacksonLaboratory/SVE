
export MKFILE_DIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SRC=src
TARGET_BIN=bin
CONFIG=$(TARGET_BIN)/tools.py
PROGRAM=$(SVE_DIR)/scripts/sve

SVE_DIR=$(shell pwd)
BROAD_BUNDLE = "ftp://ftp.broadinstitute.org/pub/svtoolkit"
GENOME_STRIP_HG38_META = $(BROAD_BUNDLE)/hg38_pre_release/Homo_sapiens_assembly38_beta1_18Aug2016.tar.gz

TARS = jre1.8.0_51 picard-tools-2.5.0 svtoolkit_2.00.1736 CNVnator_v0.3.3
SUBDIRS = bwa speedseq htslib samtools bcftools delly lumpy-sv

GCCVERSION=$(shell gcc --version | grep ^gcc | sed 's/^.* //g' | awk -F'.' '{print $1"."$2}')

# all
all: unzip_src configure build CNVnator_v0.3.3
	@test -d $(SVE_DIR)/data || tar -zxvf data.tar.gz # unzip data
	@test -d $(SVE_DIR)/$(TARGET_BIN) || mkdir $(SVE_DIR)/$(TARGET_BIN)
	$(MAKE) tools_py
	@cp $(PROGRAM) $(SVE_DIR)/$(TARGET_BIN)
.PHONY: all

unzip_src:
	@cd $(SVE_DIR)/$(SRC); \
	for module in $(TARS); do \
		if [ ! -d $$module ]; then \
			echo "- Unzip $$module"; \
			tar -zxvf $$module.tar.gz; \
		fi; \
	done; \
	cd $(SVE_DIR)

configure:
	@cd $(SVE_DIR)/$(SRC)/htslib && autoheader && autoconf && ./configure --disable-lzma && cd $(SVE_DIR)
	$(MAKE) --no-print-directory -C $(SVE_DIR)/$(SRC)/lumpy-sv

build:
	@for dir in $(SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $(SRC)/$$dir;\
	done

CNVnator_v0.3.3:
	$(MAKE) --no-print-directory -C $(SRC)/CNVnator_v0.3.3/src/samtools
	$(MAKE) --no-print-directory -C $(SRC)/CNVnator_v0.3.3/src

tools_py:
	@echo "TOOLS={}" > $(CONFIG)
	@echo "TOOLS ['SVE_HOME'] = '$(SVE_DIR)'" >> $(CONFIG)
	@echo "TOOLS ['BWA'] = '$(SVE_DIR)/$(SRC)/bwa/bwa'" >> $(CONFIG)
	@echo "TOOLS ['BWA-POSTALT'] = '$(SVE_DIR)/$(SRC)/k8 $(SVE_DIR)/$(SRC)/bwa/bwakit/bwa-postalt.js'" >> $(CONFIG)
	@echo "TOOLS ['SPEEDSEQ'] = '$(SVE_DIR)/$(SRC)/speedseq/bin/speedseq'" >> $(CONFIG)
	@echo "TOOLS ['SAMBAMBA'] = '$(SVE_DIR)/$(SRC)/speedseq/bin/sambamba'" >> $(CONFIG)
	@echo "TOOLS ['SAMTOOLS'] = '$(SVE_DIR)/$(SRC)/samtools/samtools'" >> $(CONFIG)
	@echo "TOOLS ['TABIX'] = '$(SVE_DIR)/$(SRC)/htslib/tabix'" >> $(CONFIG)
	@echo "TOOLS ['JAVA-1.8'] = '$(SVE_DIR)/$(SRC)/jre1.8.0_51/bin/java'" >> $(CONFIG)
	@echo "TOOLS ['PICARD'] = '$(SVE_DIR)/$(SRC)/picard-tools-2.5.0/picard.jar'" >> $(CONFIG)
	@echo "TOOLS ['CNVNATOR'] = '$(SVE_DIR)/$(SRC)/CNVnator_v0.3.3/src/cnvnator'" >> $(CONFIG)
	@echo "TOOLS ['CNVNATOR2VCF'] = '$(SVE_DIR)/$(SRC)/CNVnator_v0.3.3/cnvnator2VCF.pl'" >> $(CONFIG)
	@echo "TOOLS ['DELLY'] = '$(SVE_DIR)/$(SRC)/delly/src/delly'" >> $(CONFIG)
	@echo "TOOLS ['BCFTOOLS'] = '$(SVE_DIR)/$(SRC)/bcftools/bcftools'" >> $(CONFIG)
	@echo "TOOLS ['LUMPY-EXPRESS'] = '$(SVE_DIR)/$(SRC)/lumpy-sv/bin/lumpyexpress'" >> $(CONFIG)
	@echo "TOOLS ['BWA_PATH'] = '$(SVE_DIR)/$(SRC)/bwa'" >> $(CONFIG)
	@echo "TOOLS ['SAMTOOLS_PATH'] = '$(SVE_DIR)/$(SRC)/samtools'" >> $(CONFIG)
	@echo "TOOLS ['BCFTOOLS_PATH'] = '$(SVE_DIR)/$(SRC)/bcftools'" >> $(CONFIG)
	@echo "TOOLS ['HTSLIB_PATH'] = '$(SVE_DIR)/$(SRC)/htslib'" >> $(CONFIG)
	@echo "TOOLS ['JAVA-1.8_PATH'] = '$(SVE_DIR)/$(SRC)/jre1.8.0_51/bin'" >> $(CONFIG)
	@echo "TOOLS ['GENOME_STRIP_PATH'] = '$(SVE_DIR)/$(SRC)/svtoolkit_2.00.1736'" >> $(CONFIG)
	@echo "" >> $(CONFIG)
	@echo "FILES={}" >> $(CONFIG)
	@echo "FILES ['GRCH38-EXTRA'] = '$(SVE_DIR)/data/bwakit-GRCh38/hs38DH-extra.fa'"  >> $(CONFIG)
	@echo "FILES ['GRCH38-ALT'] = '$(SVE_DIR)/data/bwakit-GRCh38/hs38DH.fa.alt'"  >> $(CONFIG)
	@echo "FILES ['DELLY-HG19'] = '$(SVE_DIR)/$(SRC)/delly/excludeTemplates/human.hg19.excl.tsv'"  >> $(CONFIG)
	@echo "FILES ['DELLY-HG38'] = '$(SVE_DIR)/$(SRC)/delly/excludeTemplates/human.hg38.excl.tsv'"  >> $(CONFIG)

clean:
	@echo "- Clean up"
	@for dir in $(SUBDIRS); do \
		if [ -d $$dir ]; then \
			$(MAKE) -C $(SRC)/$$dir clean; \
		fi; \
	done
	@for module in $(TARS); do \
		if [ -d $(SRC)/$$module ]; then \
			rm -rf $(SRC)/$$module; \
		fi; \
	done
	@if [ -d $(SVE_DIR)/data ]; then \
		chmod -R 777 $(SVE_DIR)/data; \
		rm -rf $(SVE_DIR)/data; \
	fi
	@rm -rf $(TARGET_BIN)

.PHONY: clean

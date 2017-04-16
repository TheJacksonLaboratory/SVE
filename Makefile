
export MKFILE_DIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SRC=src
TARGET_BIN=bin
CONFIG=$(TARGET_BIN)/tools.py

SVE_DIR=$(shell pwd)

TARS = samtools-1.3 jre1.8.0_51 picard-tools-2.5.0
SUBDIRS = bwa speedseq samtools-1.3

# all
all:
	$(MAKE) unzip
	$(MAKE) build
	$(MAKE) config
.PHONY: all

unzip:
	@cd $(SVE_DIR)/$(SRC); \
	for module in $(TARS); do \
		if [ ! -d $$module ]; then \
			echo "- Unzip $$module"; \
			tar -zxvf $$module.tar.gz; \
		fi; \
	done; \
	cd $(SVE_DIR)

build:
	@for dir in $(SUBDIRS); do \
		echo "- Building in $$dir"; \
		$(MAKE) --no-print-directory -C $(SRC)/$$dir;\
	done
config:
	@echo "TOOLS={}" > $(CONFIG)
	@echo "TOOLS ['BWA'] = '$(SVE_DIR)/$(SRC)/bwa/bwa'" >> $(CONFIG)
	@echo "TOOLS ['BWA-POSTALT'] = '$(SVE_DIR)/$(SRC)/k8 $(SVE_DIR)/$(SRC)/bwa/bwakit/bwa-postalt.js'" >> $(CONFIG)
	@echo "TOOLS ['SPEEDSEQ'] = '$(SVE_DIR)/$(SRC)/speedseq/bin/speedseq'" >> $(CONFIG)
	@echo "TOOLS ['SAMBAMBA'] = '$(SVE_DIR)/$(SRC)/speedseq/bin/sambamba'" >> $(CONFIG)
	@echo "TOOLS ['SAMTOOLS-1.3'] = '$(SVE_DIR)/$(SRC)/samtools-1.3/samtools'" >> $(CONFIG)
	@echo "TOOLS ['JAVA-1.8'] = '$(SVE_DIR)/$(SRC)/jre1.8.0_51/bin/java'" >> $(CONFIG)
	@echo "TOOLS ['PICARD'] = '$(SVE_DIR)/$(SRC)/picard-tools-2.5.0/picard.jar'" >> $(CONFIG)
	@echo "FILES={}" >> $(CONFIG)
	@echo "FILES ['GRCH38-EXTRA'] = '$(SVE_DIR)/data/bwakit-GRCh38/hs38DH-extra.fa'"  >> $(CONFIG)
	@echo "FILES ['GRCH38-ALT'] = '$(SVE_DIR)/data/bwakit-GRCh38/hs38DH.fa.alt'"  >> $(CONFIG)

clean:
	@echo "- Clean up"
	@for dir in $(SUBDIRS); do \
		if [ -d $$dir ]; then \
			$(MAKE) --no-print-directory -C $(SRC)/$$dir clean; \
		fi; \
	done
	@for module in $(TARS); do \
		rm -rf $(SRC)/$$module; \
	done
	@rm -f $(CONFIG)

.PHONY: clean

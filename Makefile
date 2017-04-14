
export MKFILE_DIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SRC=src

SVE_DIR=$(shell pwd)

TARS = samtools-1.3 jre1.8.0_51 picard-tools-2.5.0
SUBDIRS = bwa speedseq samtools-1.3

# all
all:
	$(MAKE) unzip
	$(MAKE) build
.PHONY: all

unzip:
	@cd $(SVE_DIR)/$(SRC); \
	for module in $(TARS); do \
		if [ -f $$module.tar.gz ]; then \
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

# modules

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

.PHONY: clean

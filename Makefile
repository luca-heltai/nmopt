BUILD_DIR=jupyterbook/_build/html
DEALII_DIR=codes/dealii
DEALII_DOXYGEN_HTML=$(DEALII_DIR)/doc/html
DEALII_BOOK_HTML=$(BUILD_DIR)/doxygen/dealii

download-dealii-tag:
	@mkdir -p $(DEALII_DIR)
	curl -fsSL https://dealii.org/developer/doxygen/deal.tag -o $(DEALII_DIR)/deal.tag

build-dealii-doxygen: download-dealii-tag
	cd $(DEALII_DIR) && doxygen doc/Doxyfile

sync-doxygen-assets:
	@mkdir -p $(DEALII_BOOK_HTML)
	@rm -rf $(DEALII_BOOK_HTML)
	@mkdir -p $(DEALII_BOOK_HTML)
	cp -R $(DEALII_DOXYGEN_HTML)/. $(DEALII_BOOK_HTML)/

sync-notebook-assets:
	@mkdir -p $(BUILD_DIR)
	@if [ -d $(BUILD_DIR)/build ]; then \
		find $(BUILD_DIR)/build -maxdepth 1 -type f \( -name '*.png' -o -name '*.jpg' -o -name '*.jpeg' -o -name '*.svg' -o -name '*.webp' \) -exec cp -f {} $(BUILD_DIR)/ \; ; \
	fi

build: build-dealii-doxygen
	. ./start.sh && cd jupyterbook && jupyter book build --html --execute
	@$(MAKE) sync-doxygen-assets
	@$(MAKE) sync-notebook-assets

clean:
	. ./start.sh && cd jupyterbook && jupyter book clean -y
	rm -rf jupyterbook/_build
	find jupyterbook -type d \( -name '__pycache__' -o -name '.ipynb_checkpoints' \) -prune -exec rm -rf {} +

serve:
	make build && \
	make link-slides && \
	. ./start.sh && cd $(BUILD_DIR) && python -m http.server 3001

start: serve

dev:
	. ./start.sh && cd jupyterbook && jupyter book start

link-slides:
	echo Linking slides
	rm -rf `pwd`/$(BUILD_DIR)/slideshow
	ln -s `pwd`/jupyterbook/slides `pwd`/$(BUILD_DIR)/slideshow

copy-slides:
	echo Copying slides
	rm -rf `pwd`/$(BUILD_DIR)/slideshow
	cp -r `pwd`/jupyterbook/slides `pwd`/$(BUILD_DIR)/slideshow

slides:
	./scripts/generate_slides.sh -v

serve-slides:
	@. ./start.sh && \
	PORT=3002; \
	while lsof -iTCP:$$PORT -sTCP:LISTEN >/dev/null 2>&1; do \
		PORT=$$((PORT+1)); \
	done; \
	echo "Serving slides on http://localhost:$$PORT"; \
	cd jupyterbook/slides && python -m http.server $$PORT

all: build slides copy-slides

.PHONY: download-dealii-tag build-dealii-doxygen sync-doxygen-assets sync-notebook-assets build clean serve start dev link-slides copy-slides slides serve-slides all

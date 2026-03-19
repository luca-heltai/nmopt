BUILD_DIR=jupyterbook/_build/html

sync-notebook-assets:
	@mkdir -p $(BUILD_DIR)
	@if [ -d $(BUILD_DIR)/build ]; then \
		find $(BUILD_DIR)/build -maxdepth 1 -type f \( -name '*.png' -o -name '*.jpg' -o -name '*.jpeg' -o -name '*.svg' -o -name '*.webp' \) -exec cp -f {} $(BUILD_DIR)/ \; ; \
	fi

build:
	. ./start.sh && cd jupyterbook && jupyter book build --html --execute
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

.PHONY: sync-notebook-assets build clean serve start dev link-slides copy-slides slides serve-slides all

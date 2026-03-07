SITE_DIR=jupyterbook/_build/site/public
HTML_DIR=jupyterbook/_build/html

sync-notebook-assets:
	@mkdir -p $(HTML_DIR) $(SITE_DIR)
	@if [ -d $(HTML_DIR)/build ]; then \
		find $(HTML_DIR)/build -maxdepth 1 -type f \( -name '*.png' -o -name '*.jpg' -o -name '*.jpeg' -o -name '*.svg' -o -name '*.webp' \) -exec cp -f {} $(HTML_DIR)/ \; ; \
	fi
	@if [ -d $(SITE_DIR)/build ]; then \
		find $(SITE_DIR)/build -maxdepth 1 -type f \( -name '*.png' -o -name '*.jpg' -o -name '*.jpeg' -o -name '*.svg' -o -name '*.webp' \) -exec cp -f {} $(SITE_DIR)/ \; ; \
	fi

build:
	. ./start.sh && cd jupyterbook && jupyter book build --html
	@$(MAKE) sync-notebook-assets

clean:
	. ./start.sh && cd jupyterbook && jupyter book clean -y

serve:
	make build && \
	make link && \
	. ./start.sh && cd $(SITE_DIR) && python -m http.server 3001

start: serve

dev:
	. ./start.sh && cd jupyterbook && jupyter book start

link: 
	echo Linking slides
	rm -rf `pwd`/$(SITE_DIR)/slideshow
	ln -s `pwd`/jupyterbook/slides `pwd`/$(SITE_DIR)/slideshow

copy: 
	echo Copying slides
	rm -f `pwd`/$(SITE_DIR)/slideshow
	cp -r `pwd`/jupyterbook/slides `pwd`/$(SITE_DIR)/slideshow

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

all: build slides copy

.PHONY: sync-notebook-assets build clean serve start dev link copy slides serve-slides all

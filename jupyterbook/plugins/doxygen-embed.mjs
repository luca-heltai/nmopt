function escapeAttr(value) {
  return String(value)
    .replace(/&/g, "&amp;")
    .replace(/"/g, "&quot;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;");
}

const doxygenIframeDirective = {
  name: "doxygen-iframe",
  options: {
    src: {
      type: String,
      required: true,
    },
    title: {
      type: String,
    },
    height: {
      type: String,
    },
    width: {
      type: String,
    },
    link_text: {
      type: String,
    },
  },
  run(data) {
    const src = escapeAttr(data.options?.src ?? "");
    const title = escapeAttr(
      data.options?.title ?? "Doxygen API reference"
    );
    const height = escapeAttr(data.options?.height ?? "900px");
    const width = escapeAttr(data.options?.width ?? "100%");
    const linkText = escapeAttr(
      data.options?.link_text ?? "Open the API reference in a separate tab"
    );

    return [
      {
        type: "html",
        value: [
          '<div class="nmopt-doxygen-frame">',
          `  <iframe src="${src}" title="${title}" width="${width}" height="${height}" style="border: 1px solid #d0d7de; border-radius: 8px; background: white;"></iframe>`,
          `  <p><a href="${src}" target="_blank" rel="noopener noreferrer">${linkText}</a></p>`,
          "</div>",
        ].join("\n"),
      },
    ];
  },
};

export default {
  name: "nmopt-doxygen-embed",
  directives: [doxygenIframeDirective],
};

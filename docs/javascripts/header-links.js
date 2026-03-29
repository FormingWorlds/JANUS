function wire() {
  const homepage = "https://proteus-framework.org/";

  const logo = document.querySelector(".md-header__button.md-logo");
  if (logo) logo.href = homepage;

  const title = document.querySelector(".md-header__title[data-md-component='header-title']");
  if (title && !title.dataset.titleWired) {
    title.dataset.titleWired = "1";
    title.style.cursor = "pointer";

    // derive docs home from the first path segment (e.g. "/JANUS/" or "/Zalmoxis/"), falling back to "/" (mkdocs serve)
    const pathname = location.pathname || "/";
    const pathSegments = pathname.split("/").filter(Boolean);
    const basePath = pathSegments.length ? `/${pathSegments[0]}/` : "/";
    const docsHome = location.origin + basePath;

    title.addEventListener("click", (e) => {
      if (e.target.closest("a, button, input, label")) return;
      window.location.assign(docsHome);
    }, true);
  }
}

document.addEventListener("DOMContentLoaded", wire);
if (window.document$?.subscribe) window.document$.subscribe(wire);

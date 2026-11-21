// Sections marked .reveal fade in as they scroll into view and then stay.
// The page shows everything without JavaScript, with reduced motion, or without
// IntersectionObserver; the html.js-reveal flag is set in <head> so hidden
// sections never flash before this script runs.
(() => {
  const root = document.documentElement;
  const sections = [...document.querySelectorAll('.reveal')];
  const motion = matchMedia('(prefers-reduced-motion: reduce)');
  function showAll() { sections.forEach(section => section.classList.add('is-revealed')); }
  if (!('IntersectionObserver' in window) || motion.matches) { showAll(); return; }
  // Revealing a section also reveals every section above it, so a fast scroll or
  // a jump never leaves a hidden section behind the reader.
  function revealUpTo(section) {
    sections.slice(0, sections.indexOf(section) + 1).forEach(earlier => {
      earlier.classList.add('is-revealed');
      observer.unobserve(earlier);
    });
  }
  const observer = new IntersectionObserver(entries => entries.forEach(entry => {
    if (entry.isIntersecting || entry.boundingClientRect.bottom < 0) revealUpTo(entry.target);
  }), { rootMargin: '0px 0px -20% 0px', threshold: 0 });
  sections.forEach(section => observer.observe(section));
  // A jump to an anchor (e.g. #footprint) reveals the target at once.
  function revealHash() {
    const target = location.hash && document.getElementById(location.hash.slice(1));
    const section = target && target.closest('.reveal');
    if (section) revealUpTo(section);
  }
  addEventListener('hashchange', revealHash);
  revealHash();
  motion.addEventListener('change', () => { if (motion.matches) showAll(); });
  root.classList.add('reveal-ready');
})();

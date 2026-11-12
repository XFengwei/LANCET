(() => {
  document.querySelectorAll('.wavelength-gallery').forEach(gallery => {
  const frames = [...gallery.querySelectorAll('.gallery-frame')];
  const buttons = [...gallery.querySelectorAll('button')];
  const motion = matchMedia('(prefers-reduced-motion: reduce)');
  // data-hold (ms): a selected frame stays this long, then the five-second cycle resumes.
  // Without it, hovering or focusing the gallery holds the current frame.
  const hold = Number(gallery.dataset.hold) || 0;
  let current = 0, timer, delay;
  function show(index) {
    current = index;
    frames.forEach((frame, i) => frame.classList.toggle('is-visible', i === index));
    buttons.forEach((button, i) => button.setAttribute('aria-pressed', String(i === index)));
  }
  function interacting() {
    return !hold && (gallery.contains(document.activeElement) || gallery.matches(':hover'));
  }
  function restart() {
    clearInterval(timer);
    clearTimeout(delay);
    if (!motion.matches) timer = setInterval(() => {
      if (!document.hidden && !interacting()) show((current + 1) % frames.length);
    }, 5000);
  }
  function holdThenResume() {
    clearInterval(timer);
    clearTimeout(delay);
    if (motion.matches) return;
    delay = setTimeout(() => { show((current + 1) % frames.length); restart(); }, hold);
  }
  buttons.forEach((button, i) => button.addEventListener('click', () => {
    show(i);
    if (hold) holdThenResume(); else restart();
  }));
  motion.addEventListener('change', restart);
  Promise.all(frames.map(frame => frame.decode())).then(restart).catch(() => {});
  });
})();

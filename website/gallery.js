(() => {
  document.querySelectorAll('.wavelength-gallery, .telescope-gallery').forEach(gallery => {
  // Frames are images (.gallery-frame) or rows of link cards (.gallery-row).
  const frames = [...gallery.querySelectorAll('.gallery-frame, .gallery-row')];
  // Buttons either pick a frame (one per frame) or step by data-step (‹ = -1, › = +1).
  const steps = [...gallery.querySelectorAll('button[data-step]')];
  const buttons = [...gallery.querySelectorAll('button:not([data-step])')];
  const motion = matchMedia('(prefers-reduced-motion: reduce)');
  // data-hold (ms): a selected frame stays this long, then the five-second cycle resumes.
  // Without it, hovering or focusing the gallery holds the current frame.
  const hold = Number(gallery.dataset.hold) || 0;
  // Frames that hold links always pause while the pointer or focus is on them,
  // so a card never changes under the cursor.
  const linkArea = frames.some(frame => frame.querySelector('a')) ? frames[0].parentElement : null;
  let current = 0, timer, delay;
  function show(index) {
    current = index;
    frames.forEach((frame, i) => frame.classList.toggle('is-visible', i === index));
    buttons.forEach((button, i) => button.setAttribute('aria-pressed', String(i === index)));
  }
  function interacting() {
    if (linkArea && (linkArea.matches(':hover') || linkArea.contains(document.activeElement))) return true;
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
    delay = setTimeout(() => { restart(); if (!interacting()) show((current + 1) % frames.length); }, hold);
  }
  function choose(index) {
    show((index + frames.length) % frames.length);
    if (hold) holdThenResume(); else restart();
  }
  buttons.forEach((button, i) => button.addEventListener('click', () => choose(i)));
  steps.forEach(button => button.addEventListener('click', () => choose(current + Number(button.dataset.step))));
  motion.addEventListener('change', restart);
  const images = frames.filter(frame => frame.tagName === 'IMG');
  if (images.length) Promise.all(images.map(image => image.decode())).then(restart).catch(() => {});
  else restart();
  });
})();

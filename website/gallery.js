(() => {
  document.querySelectorAll('.wavelength-gallery').forEach(gallery => {
  const frames = [...gallery.querySelectorAll('.gallery-frame')];
  const buttons = [...gallery.querySelectorAll('button')];
  const motion = matchMedia('(prefers-reduced-motion: reduce)');
  let current = 0, timer;
  function show(index) {
    current = index;
    frames.forEach((frame, i) => frame.classList.toggle('is-visible', i === index));
    buttons.forEach((button, i) => button.setAttribute('aria-pressed', String(i === index)));
  }
  function restart() {
    clearInterval(timer);
    if (!motion.matches) timer = setInterval(() => {
      if (!document.hidden && !gallery.contains(document.activeElement) && !gallery.matches(':hover')) show((current + 1) % frames.length);
    }, 5000);
  }
  buttons.forEach((button, i) => button.addEventListener('click', () => { show(i); restart(); }));
  motion.addEventListener('change', restart);
  Promise.all(frames.map(frame => frame.decode())).then(restart).catch(() => {});
  });
})();

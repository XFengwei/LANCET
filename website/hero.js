// Every frame uses the same sky grid, crop, and parent transform.
(() => {
  const sky = document.querySelector('.hero-sky');
  const frames = [...sky.querySelectorAll('img')];
  const label = document.querySelector('#wavelength-label');
  const motion = matchMedia('(prefers-reduced-motion: reduce)');
  const labels = ['VVV · J, H, Ks', 'Spitzer + Herschel · 3.6, 5.8, 70 μm', 'Herschel · 250 μm', 'MeerKAT · 1.38 GHz'];
  let paused = motion.matches, elapsed = 0, previous = null;
  const ease = t => (1 - Math.cos(Math.PI * t)) / 2;
  function render() {
    const phase = (elapsed / 1000) % 48;
    const band = Math.floor(phase / 12);
    const local = phase % 12;
    const zoomingIn = band % 2 === 0;
    const progress = ease(Math.min(local / 10, 1));
    const scale = 1 + .3 * (zoomingIn ? progress : 1 - progress);
    sky.style.transform = `scale(${scale})`;
    const fade = ease(Math.max(0, (local - 10) / 2));
    frames.forEach((frame, i) => {
      frame.style.opacity = i === band ? 1 - fade : i === (band + 1) % 4 ? fade : 0;
    });
    label.textContent = `G316.8 / ${labels[fade < .5 ? band : (band + 1) % 4]}`;
  }
  motion.addEventListener('change', () => {
    paused = motion.matches;
    if (paused) { elapsed = 0; render(); }
  });
  function tick(now) {
    if (previous !== null && !paused && !document.hidden) elapsed += Math.min(now - previous, 100);
    previous = now;
    render();
    requestAnimationFrame(tick);
  }
  // Load all bands before advancing to prevent blank crossfades.
  Promise.all(frames.map(frame => frame.decode())).then(() => requestAnimationFrame(tick)).catch(() => {
    paused = true;
  });
})();

(() => {
  const frame = document.querySelector('.logo-story iframe');
  const observer = new IntersectionObserver(entries => {
    if (!entries.some(entry => entry.isIntersecting)) return;
    frame.src = frame.dataset.src;
    observer.disconnect();
  }, { threshold: .2 });
  observer.observe(frame);
})();

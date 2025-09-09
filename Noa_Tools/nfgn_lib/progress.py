"""Progress bar and pinned output helpers."""
import sys
import time

class ProgressBar:
    """A simple, elegant, TQDM-styled progress bar."""
    def __init__(self, total: int, prefix: str = '', width: int = 40, file=sys.stderr):
        self.total = total
        self.prefix = prefix
        self.width = width
        self.file = file
        self.start_time = time.time()
        self.done = 0
        self._draw()

    def update(self, amount: int = 1):
        self.done += amount
        self.done = min(self.done, self.total)
        self._draw()

    def set_prefix(self, prefix: str):
        self.prefix = prefix
        self._draw()

    def _draw(self):
        ratio = self.done / self.total if self.total > 0 else 0
        filled = int(self.width * ratio)
        bar = '#' * filled + '-' * (self.width - filled)
        percent = int(ratio * 100)
        
        elapsed = time.time() - self.start_time
        speed = self.done / elapsed if elapsed > 0 else 0
        
        if self.done == self.total:
            eta_str = f'{elapsed:.2f}s'
        elif speed > 0:
            eta = (self.total - self.done) / speed
            eta_str = f'{eta:.1f}s'
        else:
            eta_str = '?'
            
        self.file.write(f'{self.prefix} |{bar}| {self.done}/{self.total} [{percent}%] [{speed:.1f}it/s, ETA: {eta_str}]')
        self.file.flush()

    def finish(self):
        self.done = self.total
        self._draw()
        self.file.write('\n')
        self.file.flush()

def clear_progress_line(file=sys.stderr):
    file.write('\r' + ' ' * 120 + '\r')
    file.flush()

_active_bar = None

def progress_print_line(line: str, file=sys.stdout):
    if _active_bar:
        clear_progress_line()
    file.write(line + '\n')
    if _active_bar:
        _active_bar._draw()
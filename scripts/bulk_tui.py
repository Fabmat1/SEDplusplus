#!/usr/bin/env python3
"""bulk_tui.py -- live terminal dashboard for bulk_fit.py.

Renders a rich-based status panel (phase, batch progress + ETA, fits/s,
stars fit / remaining, CPU + RAM load) that is updated from the bulk
driver's polling loop. Falls back to plain periodic status lines when
`rich` is unavailable, stdout is not a TTY (log files, batch systems), or
--no-tui was given.

Usage (see bulk_fit.py):
    with Dashboard(enabled=...) as dash:
        dash.set_total_stars(n)
        dash.set_phase("MS-first fits", batch_total=len(cfgs))
        dash.batch_advance(n_fits, n_stars)   # from the poll loop
        dash.tick()                           # refresh rates + CPU/RAM
        dash.log("message")                   # scrolls above the panel
"""
import sys
import time
from collections import deque

try:
    import psutil
except ImportError:
    psutil = None

try:
    from rich.console import Console, Group
    from rich.live import Live
    from rich.panel import Panel
    from rich.progress import (BarColumn, Progress, ProgressColumn,
                               TaskProgressColumn, TextColumn,
                               TimeElapsedColumn, TimeRemainingColumn)
    from rich.table import Table
    from rich.text import Text
    _HAVE_RICH = True

    class _CountColumn(ProgressColumn):
        """completed/total, with '?' while the total is unknown."""
        def render(self, task):
            total = "?" if task.total is None else f"{task.total:.0f}"
            return Text(f"{task.completed:.0f}/{total}", style="cyan")
except ImportError:
    _HAVE_RICH = False

_RATE_WINDOW = 120.0        # s of history for the fits/s estimate
_PLAIN_INTERVAL = 20.0      # s between status lines in plain mode
_RENDER_INTERVAL = 0.25     # s between full panel re-renders in rich mode


def _fmt_eta(seconds):
    if seconds is None or not (seconds >= 0) or seconds == float("inf"):
        return "--:--"
    seconds = int(seconds)
    h, m, s = seconds // 3600, (seconds // 60) % 60, seconds % 60
    return f"{h}:{m:02d}:{s:02d}" if h else f"{m}:{s:02d}"


class Dashboard:
    """Campaign-level fitting dashboard. All updates are driven by the
    caller (batch_advance / tick); nothing runs in the background."""

    def __init__(self, enabled=True):
        self.total_stars = 0
        self.stars_done = 0
        self.fits_total = 0        # fits completed over the whole campaign
        self.phase = "starting"
        self.batch_total = 0
        self.batch_done = 0
        # whether batch completion of a star's configs finalizes the star
        # (False during MS-first phase 1: acceptance is only known later)
        self.count_stars = True
        self.fits_phase = True    # batch units are fits (vs prep/collect)
        self.t_start = time.time()
        self._hist = deque()       # (t, fits_total) for the rate window
        self._last_plain = 0.0
        self._last_render = 0.0
        if psutil:
            psutil.cpu_percent(None)   # prime; first call always returns 0
        self._rich = _HAVE_RICH and enabled and sys.stdout.isatty()
        self._plain = enabled and not self._rich
        self._live = None
        if self._rich:
            self._console = Console()
            self._batch_prog = Progress(
                TextColumn("[bold blue]{task.description}"),
                BarColumn(bar_width=None),
                TaskProgressColumn(),
                _CountColumn(),
                TextColumn("•"), TimeElapsedColumn(),
                TextColumn("• ETA"), TimeRemainingColumn(),
                console=self._console)
            self._batch_task = self._batch_prog.add_task("", total=1)
            self._sys_prog = Progress(
                TextColumn("{task.description}"),
                BarColumn(bar_width=24),
                TextColumn("{task.fields[label]}"),
                console=self._console)
            self._cpu_task = self._sys_prog.add_task(
                "[magenta]CPU", total=100, label="")
            self._ram_task = self._sys_prog.add_task(
                "[magenta]RAM", total=100, label="")

    # -- lifecycle ---------------------------------------------------------

    def __enter__(self):
        if self._rich:
            self._live = Live(self._render(), console=self._console,
                              refresh_per_second=4)
            self._live.__enter__()
        return self

    def __exit__(self, *exc):
        if self._live:
            self._live.update(self._render(), refresh=True)
            self._live.__exit__(*exc)
            self._live = None
        return False

    # -- state updates -----------------------------------------------------

    def set_total_stars(self, n):
        self.total_stars = n

    def set_phase(self, name, batch_total=0, count_stars=True, fits=True):
        """Start a new phase; batch progress restarts at 0/batch_total.
        fits=False: the batch counts stars, not fits (prep/collect).
        batch_total=None: total unknown (indeterminate bar, count only)."""
        self.phase = name
        self.batch_total = batch_total
        self.count_stars = count_stars
        self.fits_phase = fits
        self.batch_done = 0
        if self._rich:
            total = None if batch_total is None else max(batch_total, 1)
            self._batch_prog.reset(self._batch_task, total=total,
                                   description=name)
        self.tick(force=True)

    def batch_advance(self, n_fits, n_stars=0):
        self.batch_done += n_fits
        self.fits_total += n_fits
        self.stars_done += n_stars

    def step(self, n=1):
        """Advance the batch progress only (non-fit phases: prep, collect)."""
        self.batch_done += n

    def add_stars(self, n):
        self.stars_done += n

    def log(self, msg):
        """A scrolling message (kept above the live panel in rich mode)."""
        if self._rich:
            self._console.log(msg)
        else:
            print(msg, flush=True)

    # -- rendering ---------------------------------------------------------

    def _rate(self):
        """(windowed fits/s, overall fits/s)."""
        now = time.time()
        self._hist.append((now, self.fits_total))
        while len(self._hist) > 2 and now - self._hist[0][0] > _RATE_WINDOW:
            self._hist.popleft()
        t0, n0 = self._hist[0]
        win = (self.fits_total - n0) / (now - t0) if now > t0 else 0.0
        overall = self.fits_total / max(now - self.t_start, 1e-9)
        return win, overall

    def _sysload(self):
        if not psutil:
            return None
        vm = psutil.virtual_memory()
        return (psutil.cpu_percent(None), vm.percent,
                (vm.total - vm.available) / 2**30, vm.total / 2**30)

    def _render(self):
        win, overall = self._rate()
        stats = Table.grid(padding=(0, 2))
        stats.add_row(
            Text.assemble(("fits/s ", "bold"), (f"{win:.2f}", "green"),
                          (f"  (avg {overall:.2f})", "dim")),
            Text.assemble(("fits done ", "bold"), (f"{self.fits_total}",
                                                   "cyan")),
            Text.assemble(("stars fit ", "bold"),
                          (f"{self.stars_done}", "green"), ("/", "dim"),
                          (f"{self.total_stars}", "cyan"),
                          ("  remaining ", "bold"),
                          (f"{max(self.total_stars - self.stars_done, 0)}",
                           "yellow")))
        body = [self._batch_prog, stats]
        sl = self._sysload()
        if sl:
            cpu, ramp, used, tot = sl
            self._sys_prog.update(self._cpu_task, completed=cpu,
                                  label=f"[bold]{cpu:5.1f}%[/]")
            self._sys_prog.update(
                self._ram_task, completed=ramp,
                label=f"[bold]{ramp:5.1f}%[/] [dim]({used:.1f}/{tot:.0f} "
                      f"GB)[/]")
            body.append(self._sys_prog)
        title = (f"[bold]SED++ bulk fit[/] [dim]•[/] "
                 f"elapsed {_fmt_eta(time.time() - self.t_start)}")
        return Panel(Group(*body), title=title, border_style="blue")

    def tick(self, force=False):
        """Refresh rates, CPU/RAM and the display. Cheap enough to call per
        item in tight loops: the batch bar is updated every call (the Live
        auto-refresh picks it up), the full panel re-render is throttled."""
        if self._rich:
            self._batch_prog.update(self._batch_task,
                                    completed=self.batch_done,
                                    description=self.phase)
            now = time.time()
            if not force and now - self._last_render < _RENDER_INTERVAL:
                return
            self._last_render = now
            self._live.update(self._render())
        elif self._plain:
            now = time.time()
            if not force and now - self._last_plain < _PLAIN_INTERVAL:
                return
            self._last_plain = now
            win, overall = self._rate()
            total = "?" if self.batch_total is None else self.batch_total
            line = f"[{self.phase}] {self.batch_done}/{total}"
            if self.fits_phase:
                eta = ((self.batch_total - self.batch_done) / win
                       if win > 0 and self.batch_total else None)
                line += f" fits  {win:.2f} fits/s  ETA {_fmt_eta(eta)}"
            line += f"  stars {self.stars_done}/{self.total_stars}"
            sl = self._sysload()
            if sl:
                line += f"  CPU {sl[0]:.0f}%  RAM {sl[1]:.0f}%"
            print(line, flush=True)

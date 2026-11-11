"""Bounded local launcher controls; no scheduler or SW4 executable is used."""
import os
import select
from pathlib import Path
import sys
import tempfile
import time
import unittest

from receiver_template_startup import launch


class LauncherTimeout(unittest.TestCase):
    def test_normal_exit(self):
        with tempfile.TemporaryDirectory() as name:
            code, first = launch([sys.executable, '-S', '-c',
                                  'print("Time step 1", flush=True)'],
                                 Path(name), os.environ, timeout=5, grace=0.2)
            self.assertEqual(code, 0)
            self.assertIsNotNone(first)
            self.assertIn('Time step 1', (Path(name)/'sw4.log').read_text())

    def test_ignores_sigterm(self):
        self.check_timeout('import signal; signal.signal(signal.SIGTERM, signal.SIG_IGN); '
                           'print("Time step 1", flush=True); signal.pause()')

    def test_descendant_holds_output(self):
        # The child inherits output and ignores SIGTERM after the launcher exits.
        child = ('import os, signal, sys; from pathlib import Path; '
                 'signal.signal(signal.SIGTERM, signal.SIG_IGN); '
                 'Path("child.pid").write_text(str(os.getpid())); '
                 'os.write(int(sys.argv[1]), b"ready"); os.close(int(sys.argv[1])); signal.pause()')
        # Use an explicit readiness pipe rather than a timing sleep.
        parent = ('import os, subprocess, sys, signal; r,w=os.pipe(); '
                  'p=subprocess.Popen([sys.executable,"-S","-c",'+repr(child)+',str(w)],pass_fds=(w,)); '
                  'os.close(w); os.read(r,5); os.close(r); '
                  'print("Time step 1",flush=True); signal.pause()')
        self.check_timeout(parent)

    def check_timeout(self, script):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            start = time.monotonic()
            with self.assertRaises(TimeoutError):
                launch([sys.executable, '-S', '-c', script], root, os.environ,
                       timeout=2, grace=0.2)
            self.assertLess(time.monotonic()-start, 6)
            self.assertIn('Time step 1', (root/'sw4.log').read_text())
            if (root/'child.pid').exists() and hasattr(os, 'pidfd_open'):
                # Signal delivery is asynchronous; wait on the process-exit event.
                try:
                    fd = os.pidfd_open(int((root/'child.pid').read_text()))
                except ProcessLookupError:
                    return
                try:
                    self.assertTrue(select.select([fd], [], [], 2)[0],
                                    'task-owned child survived timeout cleanup')
                finally:
                    os.close(fd)



if __name__ == '__main__':
    unittest.main()

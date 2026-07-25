import logging
import subprocess
import unittest
from unittest.mock import patch

from sgtree._subprocess import run_capture, run_check


class RunCheckTests(unittest.TestCase):
    def test_success_returns_none(self):
        with patch("sgtree._subprocess.subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(
                args=["echo"], returncode=0
            )
            self.assertIsNone(run_check(["echo", "ok"]))
            mock_run.assert_called_once()
            kwargs = mock_run.call_args.kwargs
            self.assertIs(kwargs["check"], True)

    def test_reraises_called_process_error_after_logging(self):
        err = subprocess.CalledProcessError(returncode=2, cmd=["false"])
        with patch("sgtree._subprocess.subprocess.run", side_effect=err):
            with self.assertLogs("sgtree", level="ERROR") as cm:
                with self.assertRaises(subprocess.CalledProcessError):
                    run_check(["false"])
            self.assertTrue(any("false" in m and "2" in m for m in cm.output))

    def test_forwards_stdout_stderr_env_cwd(self):
        with patch("sgtree._subprocess.subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(args=["x"], returncode=0)
            env = {"A": "1"}
            run_check(["x"], stdout=subprocess.DEVNULL, env=env, cwd="/tmp")
            kwargs = mock_run.call_args.kwargs
            self.assertEqual(kwargs["stdout"], subprocess.DEVNULL)
            self.assertEqual(kwargs["env"], env)
            self.assertEqual(kwargs["cwd"], "/tmp")


class RunCaptureTests(unittest.TestCase):
    def test_returns_completed_process_without_raising_on_nonzero(self):
        with patch("sgtree._subprocess.subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(
                args=["x"], returncode=3, stdout="out", stderr="err"
            )
            with self.assertLogs("sgtree", level="ERROR"):
                result = run_capture(["x"])
            self.assertEqual(result.returncode, 3)
            self.assertEqual(result.stdout, "out")
            self.assertEqual(result.stderr, "err")

    def test_uses_text_mode_and_capture_output(self):
        with patch("sgtree._subprocess.subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(args=["x"], returncode=0)
            run_capture(["x"])
            kwargs = mock_run.call_args.kwargs
            self.assertTrue(kwargs["capture_output"])
            self.assertTrue(kwargs["text"])
            self.assertFalse(kwargs["check"])

    def test_integration_real_echo(self):
        # Real subprocess — proves the helper actually works end-to-end.
        result = run_capture(["echo", "hello"])
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout.strip(), "hello")

    def test_logs_stderr_on_non_zero_exit(self):
        with patch("sgtree._subprocess.subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(
                args=["x"], returncode=2, stdout="", stderr="boom details"
            )
            with self.assertLogs("sgtree", level="ERROR") as cm:
                run_capture(["x"])
            self.assertTrue(any("boom details" in m for m in cm.output))

    def test_forwards_env_and_cwd(self):
        with patch("sgtree._subprocess.subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(args=["x"], returncode=0)
            env = {"FOO": "bar"}
            run_capture(["x"], env=env, cwd="/tmp")
            kwargs = mock_run.call_args.kwargs
            self.assertEqual(kwargs["env"], env)
            self.assertEqual(kwargs["cwd"], "/tmp")


if __name__ == "__main__":
    unittest.main()

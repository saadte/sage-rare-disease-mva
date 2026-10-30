import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import BinaryIO

Command = Sequence[str]
if sys.platform == "win32":
    PROCESS_FLAGS = subprocess.CREATE_NO_WINDOW
else:
    PROCESS_FLAGS = 0


class CommandError(RuntimeError):
    pass


def run_pipeline(
    commands: Sequence[Command], *, log_path: Path, stdout: BinaryIO | None = None
) -> None:
    """Run a command pipeline and append stderr to the log.

    Send the final command's stdout to `stdout`, or to the log if None. Raise
    CommandError if any command exits with a nonzero status.
    """
    with log_path.open("ab") as log:
        processes: list[subprocess.Popen[bytes]] = []
        try:
            for index, command in enumerate(commands):
                is_last = index == len(commands) - 1
                process = subprocess.Popen(
                    command,
                    stdin=processes[-1].stdout if processes else subprocess.DEVNULL,
                    stdout=(stdout if stdout is not None else log)
                    if is_last
                    else subprocess.PIPE,
                    stderr=log,
                    creationflags=PROCESS_FLAGS,
                )
                if processes and processes[-1].stdout is not None:
                    processes[-1].stdout.close()
                processes.append(process)
            exit_codes = [process.wait() for process in processes]
        finally:
            for process in processes:
                if process.stdout is not None:
                    process.stdout.close()
                if process.poll() is None:
                    process.kill()
            for process in processes:
                process.wait()
    failures = [
        f"`{' '.join(command)}` exited with {exit_code}"
        for command, exit_code in zip(commands, exit_codes, strict=True)
        if exit_code != 0
    ]
    if failures:
        message = f"{'; '.join(failures)}; see {log_path}"
        raise CommandError(message)


def capture(command: Command) -> str:
    return subprocess.run(
        command, check=True, capture_output=True, text=True, creationflags=PROCESS_FLAGS
    ).stdout

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


class Command_error(RuntimeError):
    pass


def run_pipeline(
    commands: Sequence[Command], *, log_path: Path, stdout: BinaryIO | None = None
) -> None:
    """Pipe commands together and fail if any command fails."""
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
        raise Command_error(message)


def run(command: Command, *, log_path: Path) -> None:
    run_pipeline([command], log_path=log_path)


def capture(command: Command) -> str:
    return subprocess.run(
        command, check=True, capture_output=True, text=True, creationflags=PROCESS_FLAGS
    ).stdout

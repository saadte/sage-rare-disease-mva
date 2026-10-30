import json
import time
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path


class StageError(RuntimeError):
    pass


class StageRunner:
    """Run pipeline stages and reuse completed outputs.

    Skip a stage when its `.stages` marker matches the parameters and output paths and
    all outputs exist. Otherwise, delete its marker and outputs before running it.
    Write stage logs to `logs`.
    """

    def __init__(self, directory: Path) -> None:
        self._marker_directory = directory / ".stages"
        self._log_directory = directory / "logs"
        self._marker_directory.mkdir(parents=True, exist_ok=True)
        self._log_directory.mkdir(parents=True, exist_ok=True)

    def run(
        self,
        name: str,
        action: Callable[[Path], None],
        *,
        outputs: Sequence[Path],
        parameters: Mapping[str, object],
    ) -> None:
        """Run `action` with a log path unless matching outputs are already complete.

        Raise StageError when an output is missing after the action returns.
        """
        marker_path = self._marker_directory / f"{name}.json"
        marker = json.dumps(
            {
                "parameters": parameters,
                "outputs": [str(path.resolve()) for path in outputs],
            },
            sort_keys=True,
        )
        is_complete = (
            marker_path.exists()
            and marker_path.read_text(encoding="utf8") == marker
            and all(output.exists() for output in outputs)
        )
        if is_complete:
            print(f"[skip] {name}", flush=True)
            return

        print(f"[run ] {name}", flush=True)
        marker_path.unlink(missing_ok=True)
        for output in outputs:
            output.unlink(missing_ok=True)
        log_path = self._log_directory / f"{name}.log"
        log_path.write_bytes(b"")
        started = time.monotonic()
        action(log_path)
        elapsed_seconds = time.monotonic() - started
        missing = [str(output) for output in outputs if not output.exists()]
        if missing:
            message = f"stage {name} finished without writing {missing}"
            raise StageError(message)

        temporary_path = marker_path.with_suffix(".tmp")
        temporary_path.write_text(marker, encoding="utf8")
        temporary_path.replace(marker_path)
        print(f"[done] {name} in {elapsed_seconds:.0f}s", flush=True)

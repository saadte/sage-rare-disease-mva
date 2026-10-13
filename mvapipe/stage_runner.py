import json
import time
from collections.abc import Callable, Mapping, Sequence
from datetime import UTC, datetime
from pathlib import Path
from typing import TypedDict

Stage_result = dict[str, int]


class _Marker(TypedDict):
    stage: str
    parameters: str
    seconds: float
    result: Stage_result
    finished: str


class Stage_error(RuntimeError):
    pass


class Stage_runner:
    """Resume completed stages with matching parameters and outputs."""

    def __init__(self, directory: Path) -> None:
        self._marker_directory = directory / ".stages"
        self._log_directory = directory / "logs"
        self._marker_directory.mkdir(parents=True, exist_ok=True)
        self._log_directory.mkdir(parents=True, exist_ok=True)
        self.seconds: dict[str, float] = {}
        self.results: dict[str, Stage_result] = {}

    def run(
        self,
        name: str,
        action: Callable[[Path], Stage_result],
        *,
        outputs: Sequence[Path],
        parameters: Mapping[str, object],
    ) -> None:
        marker_path = self._marker_directory / f"{name}.json"
        serialised_parameters = _serialise_parameters(
            {
                "parameters": parameters,
                "outputs": [str(path.resolve()) for path in outputs],
            }
        )
        if marker_path.exists():
            marker: _Marker = json.loads(marker_path.read_text(encoding="utf8"))
            is_current = marker["parameters"] == serialised_parameters
            if is_current and all(output.exists() for output in outputs):
                print(f"[skip] {name}", flush=True)
                self.seconds[name] = marker["seconds"]
                self.results[name] = marker["result"]
                return

        print(f"[run ] {name}", flush=True)
        marker_path.unlink(missing_ok=True)
        for output in outputs:
            output.unlink(missing_ok=True)
        log_path = self._log_directory / f"{name}.log"
        log_path.write_bytes(b"")
        started = time.monotonic()
        result = action(log_path)
        seconds = round(time.monotonic() - started, 1)
        missing = [str(output) for output in outputs if not output.exists()]
        if missing:
            message = f"stage {name} finished without writing {missing}"
            raise Stage_error(message)

        marker = {
            "stage": name,
            "parameters": serialised_parameters,
            "seconds": seconds,
            "result": result,
            "finished": datetime.now(UTC).isoformat(timespec="seconds"),
        }
        temporary_path = marker_path.with_suffix(".tmp")
        temporary_path.write_text(json.dumps(marker, indent=1), encoding="utf8")
        temporary_path.replace(marker_path)
        self.seconds[name] = seconds
        self.results[name] = result
        print(f"[done] {name} in {seconds:.0f}s", flush=True)


def _serialise_parameters(parameters: Mapping[str, object]) -> str:
    return json.dumps(parameters, sort_keys=True)

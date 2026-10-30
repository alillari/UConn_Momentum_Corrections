"""Local ROOT-skim adapter with logical field mappings and stable row indices."""

from pathlib import Path

import uproot


def root_files(directory):
    files = sorted(Path(directory).resolve().glob("*.root"))
    if not files:
        raise ValueError("No ROOT input files")
    return files


def require_separate_output(directory, output):
    directory, output = Path(directory).resolve(), Path(output).resolve()
    if output == directory or directory in output.parents:
        raise ValueError("Output must be outside the input data directory")
    output.mkdir(parents=True, exist_ok=True)
    return output


def read_file(path, tree_name, mapping):
    before = path.stat()
    with uproot.open(
        path,
        handler=uproot.source.file.MemmapSource,
        decompression_executor=uproot.source.futures.TrivialExecutor(),
        interpretation_executor=uproot.source.futures.TrivialExecutor(),
    ) as f:
        tree = f[tree_name]
        missing = set(mapping.values()) - set(tree.keys())
        if missing:
            raise ValueError(f"Missing branches in {path.name}: {sorted(missing)}")
        start = 0
        for arrays in tree.iterate(
            list(mapping.values()), step_size="40 MB", library="np"
        ):
            data = {key: arrays[value] for key, value in mapping.items()}
            yield start, data
            start += len(next(iter(data.values())))
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"Input changed during processing: {path}")

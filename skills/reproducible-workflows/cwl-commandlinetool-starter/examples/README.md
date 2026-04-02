# Example

- `hello.cwl`: minimal `CommandLineTool` that writes the requested message to `greeting.txt`.
- `hello-job.yml`: sample input object with a greeting string.
- `../scripts/run_cwl_hello.py`: wrapper that validates the CWL document, runs `cwltool`, and writes a JSON summary.

# Installation

## Editable install

```bash
cd TorchGWAS
python -m pip install -e .
```

## Notes

- TorchGWAS runs on CPU by default in this environment.
- If CUDA is available, future releases can route heavy scans to GPU without changing the CLI.
- The package currently requires `numpy`, `scipy`, and `torch`.


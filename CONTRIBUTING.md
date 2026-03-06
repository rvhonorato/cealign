# Contributing

Thank you for your interest in contributing to `cealign`!

## Getting started

1. Fork the repository and clone your fork
2. Ensure you have a recent stable Rust toolchain (>= 1.87)
3. Install system dependencies (required for the `plot` feature on Linux):
   ```bash
   sudo apt-get install libfontconfig-dev
   ```

## Development workflow

Build and test:

```bash
cargo build
cargo test --all-features
```

Check formatting and lints before opening a PR:

```bash
cargo fmt --all
cargo clippy --all-features -- -D warnings
```

## Submitting changes

- Open an issue first for non-trivial changes so we can discuss the approach
- Keep PRs focused — one feature or fix per PR
- Add or update tests for any changed behaviour
- Ensure `cargo fmt`, `cargo clippy`, and `cargo test` all pass

## Reporting bugs

Please use the [bug report template](.github/ISSUE_TEMPLATE/bug_report.md) and include a minimal reproducible example where possible.

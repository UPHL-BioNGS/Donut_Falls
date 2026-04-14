# Contributing to Donut Falls

First off, thank you for considering contributing to Donut Falls! It’s people like you that make this a great tool for the public health community.

## How Can I Contribute?

### Reporting Bugs
* **Check the FAQ and Wiki:** Your issue might already be addressed there.
* **Search Existing Issues:** Before opening a new one, check if someone else has already reported it.
* **Explain what happened:** When opening an issue, provide as much detail as possible (Nextflow version, command used, and error logs).

### Suggesting Enhancements
We are always looking to add new assemblers or polishing tools! 
* Open an issue with the tag `enhancement`.
* Explain why this tool would be beneficial (e.g., "It handles high-accuracy ONT reads better than Flye").

### Pull Requests
1. **Fork the repo** and create your branch from `main`.
2. **If you've added a new tool**, ensure you've updated the `nextflow.config` and added a test case in `.github/workflows/`.
3. **Ensure tests pass.** Run `nextflow run . -profile docker,test` to make sure the core logic still works.
4. **Update documentation.** If you change a parameter, update the `nextflow_schema.json` and the Wiki.
5. **Submit the PR** with a clear description of the changes.

## Development Environment
We recommend using **Docker** or **Singularity** for development to ensure dependency versions remain consistent with the `staphb` containers used in the workflow.

## Style Guide
* **Nextflow:** Follow [nf-core](https://nf-co.re/) DSL2 best practices where possible.
* **Commits:** Use descriptive commit messages (e.g., `feat: add dragonflye assembler` instead of `stuff`).

## Questions?
Feel free to open an issue or reach out to Erin Young at eriny@utah.gov.

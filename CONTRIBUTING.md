# Contributing guidelines

## Development

### Building the documentation

The documentation is written in [markdown](https://www.markdownguide.org/basic-syntax/), and uses [mkdocs](https://www.mkdocs.org/) to generate the pages.

To build the documentation for yourself:

```console
pip install -e .[docs]
mkdocs serve
```

You can find the documentation source in the [docs](https://github.com/FormingWorlds/JANUS/tree/main/docs) directory.
If you are adding new pages, make sure to update the listing in the [`mkdocs.yml`](https://github.com/FormingWorlds/JANUS/blob/main/mkdocs.yml) under the `nav` entry.

### Running tests

JANUS uses [pytest](https://docs.pytest.org/en/latest/) to run the tests. You can run the tests for yourself using:

```console
pytest
```

To check coverage:

```console
coverage run -m pytest
coverage report  # to output to terminal
coverage html    # to generate html report
```
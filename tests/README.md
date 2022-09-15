
First of all, make sure you have installed the developer dependencies through conda.
Then activate the conda environment `small_peptides`.
To run the tests in a repository in your home directory (where you avoid errors due to resource allocation etc..), you type: `pytest --tag integration_test_<profile> --basetemp=</path/to/suitable/dir/>`.

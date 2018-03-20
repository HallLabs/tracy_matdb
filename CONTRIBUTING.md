# Contributing to matdb

Thank you for taking the time to contribute!

The following are guidlines to help contribute to the matdb python
package. These are guidlines are open to adaptation upon suggestion
and at the contributors discretion.

#### Table of Contents

[How To Contribute](#how-to-contribute)
  * [Reporting Bugs](#reporting-bugs)
  * [Suggesting Enhancements](#suggesting-enhancements)
  * [Adding Database Types](#adding-database-types)
  * [Adding Calculators](#adding-calculators)
  * [Testing](#testing)
  * [Pull Requests](#pull-requests)

[Styleguides](#styleguides)
  * [Git Commit Messages](#git-commit-messages)
  * [History Messages](#history-messages)
  * [Python Styleguide](#python-styleguide)
  * [Documentation Styleguide](#documentation-styleguide)

[Additional Notes](#additional-notes)
  * [Issue and Pull Request Labels](#issue-and-pull-request-labels)
  * [Attribution](#attribution)

## How To Contribute

### Reporting Bugs

This section guides you through submitting a bug report for
'matdb'. Following these guidlines helps developers understand your
report and reproduce the issue.

> **Note:** If you find a **Closed** issue that seems like it is the
    same thing that you're experiencing, open a new issue and include
    a link to the original issue in the body of your new one.

#### How do I Submit A (Good) Bug Report?


### Suggesting Enhancements

This section guides you through submitting an enhancement suggestion
for 'matdb', including completely new features and minor improvements to
existing functionality. 

#### How Do I Submit A (Good) Enhancement Suggestion?

Enhancement suggestions are tracked as [GitHub
issues](https://guides.github.com/features/issues/). To create a new
enhancement suggestion open an issue and provide the following
information:

* **Use a clear and descriptive title** for the issue to identify the
    suggestion.
* **Provide a step-by-step description of the suggested enhancement**
    in as many details as possible.
* **Describe the current behavior**, if applicable, and **explain
    which behavior you expected to see instead** and why.

* **Explain why this enhancement would be useful** to 'matdb' users.

### Adding Database Types

If you want to contribute a new type of database to the repository
then follow the guidelines found [here](instructions/Contribute_Database.md).

### Adding Calculators

If you want to contribute a new type of database to the repository
then follow the guidelines found
[here](instructions/Contribute_Calculator.md).

### Testing

Any new subroutine or functionality for 'matdb' must be unit tested
before a pull request will be accepted. New unit tests should:

* Be placed in the tests/(sub-package name) folder within a file titled test_'module
  name'.py where 'module name' is the module being tested.
* Any test input or output not in the test_*.py file should be placed
  in a folder within the tests/data folder.
* Unit tests must be designed to cover any new code written, i.e., the
  new unit tests must maintain %100 test coverage.

### Pull Requests

* Do not include issue numbers in the PR title
* Follow the [Python](#python-styleguide) styleguide.
* Document new code based on the [Documentation Styleguide](#documentation-styleguide).
* List relavent issue numbers in the PR body.
* PR's will only be accepted if they pass all unit tests, don't lower
  the codecoverage, quantified code finds no new issues in the PR, and
  the HISTORY.md and README.md have been updated to reflect changes.

## Styleguides

### Git Commit Messages

* Git commit messages should be short and reference the new version
  number, found in HISTORY.md, for this commit. All details of the
  commit changes should be recorded in HISTORY.md.

### History Messages

* Reference issues and pull requests liberally.
* Detail which modules were changed and how.

### Python Styleguide

'Matdb' follows the [Google Python Style](https://google.github.io/styleguide/pyguide.html).

### Documentation Styleguide

All new subroutines must be documented in order to keep the code maintainable and easy to use.

* Use [Google Style Python](http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html).

## Additional Notes

### Issue and Pull Request Labels

This section lists the labels used to help track and manage issues and pull requests.

[GitHub search](https://help.github.com/articles/searching-issues/)
makes it easy to use labels for finding groups of issues or pull
requests you're interested in. To help you find issues and pull
requests, each label is listed with search links for finding open
items with that label in `matdb` only and also across all Atom
repositories. We encourage you to read about [other search
filters](https://help.github.com/articles/searching-issues/) which
will help you write more focused queries.

Please open an issue if you have suggestions for new labels.

#### Type of Issue and Issue State

| Label name | `matdb` :mag_right: | Description |
| --- | --- | --- |
| `enhancement` | [search][search-label-enhancement]  | Feature requests. |
| `bug` | [search][search-label-bug] | Confirmed bugs or reports that are very likely to be bugs. |
| `question` | [search][search-label-question] | Questions more than bug reports or feature requests (e.g. how do I do X). |
| `feedback` | [search][search-label-feedback] | General feedback more than bug reports or feature requests. |
| `help-wanted` | [search][search-label-help-wanted] | The matdb development team would appreciate help from the community in resolving these issues. |
| `more-information-needed` | [search][search-label-more-information-needed] | More information needs to be collected about these problems or feature requests (e.g. steps to reproduce). |
| `needs-reproduction` | [search][search-label-needs-reproduction] | Likely bugs, but haven't been reliably reproduced. |
| `duplicate` | [search][search-label-duplicate] | Issues which are duplicates of other issues, i.e. they have been reported before. |
| `wontfix` | [search][search-label-wontfix] | The matdb development team has decided not to fix these issues for now, either because they're working as intended or for some other reason. |
| `invalid` | [search][search-label-invalid] | Issues which aren't valid (e.g. user errors). |

### Attribution

This document was adapted from [Contribute].

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
[Contribute]: https://github.com/atom/atom/blob/master/CONTRIBUTING.md
    
[search-label-enhancement]: https://github.com/rosenbrockc/matdb/issues?q=is%3Aissue+is%3Aopen+label%3Aenhancement
[search-label-bug]: https://github.com/rosenbrockc/matdb/issues?utf8=%E2%9C%93&q=is%3Aissue%20is%3Aopen%20label%3Abug
[search-label-question]: https://github.com/rosenbrockc/matdb/issues?utf8=%E2%9C%93&q=is%3Aissue%20is%3Aopen%20label%3Aquestion
[search-label-feedback]: https://github.com/rosenbrockc/matdb/issues?utf8=%E2%9C%93&q=is%3Aissue%20is%3Aopen%20label%3Afeedback
[search-label-help-wanted]: https://github.com/rosenbrockc/matdb/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22
[search-label-beginner]: https://github.com/rosenbrockc/matdb/issues?utf8=%E2%9C%93&q=is%3Aissue%20is%3Aopen%20label%3A%22beginner%22%20
[search-label-more-information-needed]: https://github.com/rosenbrockc/matdb/issues?q=is%3Aopen+is%3Aissue+label%3A%22more+information+needed%22
[search-label-needs-reproduction]: https://github.com/rosenbrockc/matdb/issues?q=is%3Aopen+is%3Aissue+label%3A%22needs+reproduction%22
[search-label-duplicate]: https://github.com/rosenbrockc/matdb/issues?q=is%3Aopen+is%3Aissue+label%3Aduplicate
[search-label-wontfix]: https://github.com/rosenbrockc/matdb/issues?q=is%3Aopen+is%3Aissue+label%3Awontfix
[search-label-invalid]: https://github.com/rosenbrockc/matdb/issues?q=is%3Aopen+is%3Aissue+label%3Ainvalid

[help-wanted]: https://github.com/rosenbrockc/matdb/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22
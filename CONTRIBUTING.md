## How to contribute to the catalog itself?
To contribute to the content of the catalog (addition of lenses, fixing errors, etc), 
please fork this repository and create a new branch to work on. 
Once you are satisfied with it, create a pull request [here](https://github.com/lenscat/lenscat/pulls).

We have included a simple routine in the package, as `lenscat.utils.check_possible_duplicates()`, to generate a report of possible duplicates in a given catalog. Simply run
```python
import lenscat; from lenscat.utils import check_possible_duplicates

check_possible_duplicates(lenscat.catalog)
```
with your modified catalog. In the pull request, please either show the output of this routine, or explicitly mention that you had run this routine and found no duplicates.

***Please remember to assign the label `catalog` to the pull request***. By doing so, once the pull request is merged to the main branch, 
a [github action](https://github.com/lenscat/lenscat/blob/main/.github/workflows/changes-in-catalog.yml) will be triggered, where it will perform the following actions _automatically_:
- advance the build version by one
- re-generate the visualization plot as shown in the README
- create a release for the new version
- upload the new version to PyPI

Therefore, there is no need to change the version manually when only the content of the catalog is being modified in a pull request.

## How to contribute to the `python` codebase?
To contribute to the `python` codebase (new feature, bug fix, etc), 
please fork this repository and create a new branch to work on. 
Once you are satisfied with it, create a pull request [here](https://github.com/lenscat/lenscat/pulls).

If the code change is significant, please manually advance the version (major or minor version as you see fit) in [this file](https://github.com/lenscat/lenscat/blob/main/lenscat/_version.py). After the pull request is merged to the main branch, the _approver of the pull request_ should create and publish a new release on github, and after that a [github action](https://github.com/lenscat/lenscat/blob/main/.github/workflows/python-publish.yml) will be triggered and the new version will be uploaded to PyPI automatically.

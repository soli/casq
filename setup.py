# noqa: D100
# pylint: disable=missing-docstring
from setuptools import find_packages, setup  # type: ignore

with open("README.rst", "r") as file:
    README = file.read()

with open("casq/__init__.py", "r") as file:
    for line in file.readlines():
        if line.startswith("version ="):
            version = line.split('"')[1]

setup(
    name="casq",
    version=version,
    description="CaSQ: Celldesigner as Sbml-Qual",
    long_description=README,
    long_description_content_type="text/x-rst",
    url="https://lifeware.inria.fr/~soliman/post/casq/",
    project_urls={
        "Documentation": "https://readthedocs.org/projects/casq/",
        "Code": "https://gitlab.inria.fr/soliman/casq/",
        "Mirror": "https://github.com/soli/casq",
    },
    license="GPLv3",
    packages=find_packages(),
    author="Sylvain Soliman",
    author_email="Sylvain.Soliman@inria.fr",
    install_requires=["networkx>=2.2", "loguru>=0.2.5"],
    python_requires=">=3.5",
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    extras_require={"dev": ["black", "check-manifest", "flake8-mypy"]},
    include_package_data=True,
    zip_safe=False,
    entry_points={"console_scripts": ["casq = casq.celldesigner2qual:main"]},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)

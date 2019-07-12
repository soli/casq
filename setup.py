from setuptools import find_packages, setup  # type: ignore

with open("README.md", "r") as file:
    README = file.read()

with open("casq/__init__.py", "r") as file:
    for line in file.readlines():
        if line.startswith("version ="):
            version = line.split()[2]

setup(
    name="casq",
    version=version,
    description="CaSQ: Celldesigner as Sbml-Qual",
    long_description=README,
    long_description_content="text/markdown",
    url="https://lifeware.inria.fr/~soliman/post/casq/",
    project_urls={"Source Code": "https://gitlab.inria.fr/soliman/casq/"},
    license="MIT",
    packages=find_packages(),
    author="Sylvain Soliman",
    author_email="Sylvain.Soliman@inria.fr",
    install_requires=["networkx>=2.2", "loguru>=0.2.5"],
    python_requires=">=3.5",
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    extras_require={"dev": ["black", "check-manifest", "flake8-mypy"]},
    include_package_data=True,
    data_files=[("utils", ["utils/validate.sh"])],
    zip_safe=False,
    entry_points={"console_scripts": ["casq = casq.celldesigner2qual:main"]},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

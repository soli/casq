from setuptools import find_packages, setup

with open("README.md", "r") as file:
    readme = file.read()

with open("casq/__init__.py", "r") as file:
    for line in file.readlines():
        if line.startswith("version ="):
            version = line.split()[2]

setup(
    name="casq",
    version=version,
    description="CaSQ: Celldesigner as Sbml-Qual",
    long_description=readme,
    url="https://lifeware.inria.fr/~soliman/post/casq/",
    project_urls={
        "Source Code": "https://gitlab.inria.fr/soliman/sbgnpd2sbmlq",
    },
    license="MIT",
    packages=find_packages(),
    author="Sylvain Soliman",
    author_email="Sylvain.Soliman@inria.fr",
    install_requires=["networkx>=2.2"],
    python_requires=">=3.5",
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    include_package_data=True,
    zip_safe=False,
    entry_points={
        "console_scripts": ["celldesigner2qual=casq.celldesigner2qual.main"],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="casq",
    version="0.1.0",
    description="CaSQ: Celldesigner as Sbml-Qual",
    long_description=long_description,
    url="https://lifeware.inria.fr/~soliman/post/casq/",
    license="MIT",
    packages=["casq"],
    author="Sylvain Soliman",
    author_email="Sylvain.Soliman@inria.fr",
    install_requires=["networkx"],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

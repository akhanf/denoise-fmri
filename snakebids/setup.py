import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="snakebids", 
    version="0.0.1",
    author="Ali Khan",
    author_email="alik@robarts.ca",
    description="A bids formatter for snakemake workflows",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/khanlab/snakebids",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7'
)

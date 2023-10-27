import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="optgraphstate",
    version="0.3.0",
    author="Seok-Hyung Lee",
    author_email="sh.lee1524@gmail.com",
    description="Graph-theoretical optimization of fusion-based graph state generation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/seokhyung-lee/OptGraphState",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

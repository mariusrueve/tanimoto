from setuptools import setup

def setup_package():
    """
    Set up the Tanimoto package.

    :return: None
    """
    setup(
        name="tanimoto",
        version="0.1.0",
        author="Marius Rüve",
        author_email="nope",
        description="A CLI tool for calculating Tanimoto similarity",
        long_description=open('README.md').read(),
        long_description_content_type="text/markdown",
        url="https://github.com/mariusrueve/tanimoto",
        packages=["tanimoto"],
        entry_points={
            "console_scripts": [
                "tanimoto-cli = tanimoto.cli:main"
            ]
        },
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ],
        python_requires=">=3.6",
    )

if __name__ == '__main__':
    setup_package()

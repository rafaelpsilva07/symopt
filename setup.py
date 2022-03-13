import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="symopt",
    version="1.0.0",
    install_requires=["sympy>=1.9"],
    author="Rafael Pereira",
    author_email="rafaelpsilva07@gmail.com",
    description="A support library to calculate optimization using sympy",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
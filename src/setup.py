from setuptools import setup, find_packages

setup(
    name="cca_manifolds",
    version="1.0.0",
    description="A project for analyzing gene expression patterns and drug sensitivity using CCA.",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/your-repo-url",  # Optional: Replace with your repository URL
    packages=find_packages(
        where="src"
    ),  # Automatically find packages in the src directory
    package_dir={
        "": "src"
    },  # Tell setuptools to look for packages in the src directory
    python_requires=">=3.7",  # Specify the Python version
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "scikit-learn",  # Add other dependencies your project requires
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)

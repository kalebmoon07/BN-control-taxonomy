from setuptools import setup, find_packages

setup(
    name="bntaxonomy",
    version="0.1.0",
    packages=find_packages(where="src", include=["bntaxonomy", "bntaxonomy.*"]),
    package_dir={"": "src"},
    install_requires=[],
    # entry_points={"console_scripts": []},
    author="",
    author_email="",
    description="",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="",
    classifiers=[],
    python_requires=">=3.10",
    zip_safe=False,
)

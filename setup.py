from setuptools import setup

setup(
    name="daltools",
    author="Olav Vahtras",
    author_email="olav.vahtras@gmail.com",
    version="1.0",
    url="https://github.com/vahtras/daltools",
    packages=["daltools"],
    install_requires=["util"],
    dependency_links=["https://github.com/vahtras/util.git@master#egg=util"],
    )

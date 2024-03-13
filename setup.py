import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SSWYieldModelingTutorial",
    version="0.0.1",
    author="Dmitry Savransky",
    author_email="ds264@cornell.edu",
    description="An interactive tutorial on yield modeling.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dsavransky/SSWYieldModelingTutorial",
    packages=['SSWYieldModelingTutorial'],
    include_package_data=True,
    install_requires=[
          'scipy',
          'numpy',
          'astropy',
          'sympy',
          'notebook',
          'jupyterlab',
          'matplotlib',
          'ipympl',
          'ipywidgets'
    ],
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)

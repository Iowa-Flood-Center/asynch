[![Build Status](https://travis-ci.org/Iowa-Flood-Center/asynch.svg?branch=master)](https://travis-ci.org/Iowa-Flood-Center/asynch)

# ASYNCH

A numerical library for solving differential equations with a tree structure. Emphasis is given to hillslope-link river basin models.

## Documentation

The documentation is available [here](http://asynch.readthedocs.io/). Thank you to the people running Read the Docs for such an excellent service.

The source for the documentation is in the `docs` folder. Here are the instructions to built and read it locally. The documentation is built with [Doxygen](http://www.doxygen.org/) and [Sphynx](http://www.sphinx-doc.org). The sphynx template is from [ReadtheDocs](https://docs.readthedocs.io). [Breathe](https://breathe.readthedocs.io) provides a bridge between the Sphinx and Doxygen documentation systems.

    pip install sphynx sphinx-autobuild breathe
    apt-get doxygen

    cd docs  
    doxygen api.dox
    doxygen devel.dox
    make html

The html documentation is generated in `docs/.build/html`.

## Testing

Asynch doesn't have a good test covergage at the moment but the unit test framework is in place.

import site
import sys

site.ENABLE_USER_SITE = "--user" in sys.argv[1:]


#!/usr/bin/env python
import setuptools

if __name__ == "__main__":
    setuptools.setup()

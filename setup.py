from setuptools import setup

packages_wannierpy = ['wannierpy']

if __name__ == '__main__':
    setup(name='wannierpy',
          version='0.1',
          description='Helper funtions to handle wannier input/output files with python.',
          author='Henrique Miranda',
          author_email='miranda.henrique@gmail.com',
          packages=packages_wannierpy,
          )
